#pragma once
#include <TFile.h>
#include <TTree.h>
#include "common/Logger.hpp"
#include <memory>
#include <string>
#include <unordered_map>
#include <typeindex>
#include <stdexcept>
#include <type_traits>
#include <cstddef>
#include "common/RunContext.hpp"


class RootOutput {
public:
    RootOutput(const std::string& filename,
             const std::string& treename = "events")
    : m_file(TFile::Open(filename.c_str(), "RECREATE")),
      m_tree(std::make_unique<TTree>(treename.c_str(), treename.c_str()))
    {
        if (!m_file || m_file->IsZombie()) {
            LOG_ERROR("Failed to open ROOT file: {}", filename);
            throw std::runtime_error("Failed to open ROOT file: " + filename);
        }
        // Attach tree to file so baskets can be flushed to disk during Fill()
        m_tree->SetDirectory(m_file.get());

        // Optional: reduce in-memory growth for long jobs (tune)
        m_tree->SetAutoFlush(1000);          // flush baskets every ~1000 entries
        m_tree->SetAutoSave(-50'000'000);    // autosave every ~50MB

        m_run_context_tree = std::make_unique<TTree>("run_context", "run_context");
        m_run_context_tree->SetDirectory(m_file.get());
    }

    ~RootOutput() {
        if (m_file) {
            m_file->cd();
            m_tree->Write();
            m_run_context_tree->Write();

            // Detach before closing so the file doesn't try to own/delete it.
            m_tree->SetDirectory(nullptr);
            m_run_context_tree->SetDirectory(nullptr);

            m_file->Close();
        }
    }

    TTree* tree() { return m_tree.get(); }

    /**
    * Get or create a branch buffer of type T.
    *
    * Usage:
    *   auto* buf = out.get_or_make_buffer<std::vector<float>>("RecoHits.E");
    *   buf->clear();
    *   buf->push_back(...);
    */
    template <class T>
    T* get_or_make_buffer(const std::string& branch_name) {
        auto* tree = (m_active_tree == ActiveTree::Event) ? m_tree.get() : m_run_context_tree.get();
        auto& branch_types = (m_active_tree == ActiveTree::Event) ? m_branch_types : m_run_context_branch_types;
        auto& buffers = (m_active_tree == ActiveTree::Event) ? m_buffers : m_run_context_buffers;
        const std::type_index want(typeid(T));

        auto it_type = branch_types.find(branch_name);
        if (it_type != branch_types.end() && it_type->second != want) {
        LOG_ERROR("Branch '{}' requested with different type. existing={}, requested={}",
                    branch_name, it_type->second.name(), want.name());
        throw std::runtime_error("Branch type mismatch: " + branch_name);
        }

        auto it = buffers.find(branch_name);
        if (it == buffers.end()) {
            auto holder = std::make_unique<Holder<T>>();
            T* raw_ptr = &holder->value;

            tree->Branch(branch_name.c_str(), raw_ptr);
            LOG_DEBUG("Created branch '{}' of type {}", branch_name, want.name());
            branch_types.emplace(branch_name, want);
            buffers.emplace(branch_name, std::move(holder));
            return raw_ptr;
        }

        return static_cast<Holder<T>*>(it->second.get())->ptr();
    }

    void fill() {
        m_tree->Fill();
        ++m_fill_count;
        const bool allow_shrink = (m_fill_count % m_shrink_interval == 0);
        for (auto& kv : m_buffers) {
            kv.second->reset_after_fill(allow_shrink);
        }
    }

    void write_event_counter(long long event_counter) {
        const auto active_tree = m_active_tree;
        m_active_tree = ActiveTree::Event;
        *get_or_make_buffer<long long>("EventMeta.event_counter") = event_counter;
        m_active_tree = active_tree;
    }

    void write_run_context(const RunContext& run_context, long long start_event, long long end_event) {
        const auto active_tree = m_active_tree;
        m_active_tree = ActiveTree::RunContext;
        *get_or_make_buffer<long long>("ContextRange.start_event") = start_event;
        *get_or_make_buffer<long long>("ContextRange.end_event") = end_event;

        const auto& run_cfg_desc = describe((const RunConfig*)nullptr);
        for (const auto& f : run_cfg_desc) f.write(&run_context.config, *this, "RunConfig");

        const auto& cond_desc = describe((const ConditionStore*)nullptr);
        for (const auto& f : cond_desc) f.write(&run_context.conditions, *this, "ConditionStore");

        m_run_context_tree->Fill();
        for (auto& kv : m_run_context_buffers) {
            kv.second->reset_after_fill(true);
        }
        m_active_tree = active_tree;
    }

    private:
        enum class ActiveTree {
            Event,
            RunContext
        };

        struct IHolder {
            virtual ~IHolder() = default;
            virtual void reset_after_fill(bool allow_shrink) = 0;
        };

        template <class T>
        struct Holder final : IHolder {
            T value{};
            T* ptr() { return &value; }

            void reset_after_fill(bool allow_shrink) override {
                reset_value(value, allow_shrink);
            }

        private:
            template <class X>
            static constexpr bool has_clear_v = requires(X& x) {
                x.clear();
            };

            template <class X>
            static constexpr bool has_capacity_v = requires(const X& x) {
                x.capacity();
            };

            template <class X>
            static constexpr bool has_shrink_to_fit_v = requires(X& x) {
                x.shrink_to_fit();
            };

            static void reset_value(T& v, bool allow_shrink) {
                if constexpr (has_clear_v<T>) {
                    v.clear();
                    if constexpr (has_capacity_v<T> && has_shrink_to_fit_v<T>) {
                        constexpr std::size_t max_retained_capacity = 4096;
                        if (allow_shrink && v.capacity() > max_retained_capacity) {
                            T tmp;
                            v.swap(tmp);
                        }
                    }
                }
            }
        };

        std::unique_ptr<TFile> m_file;
        std::unique_ptr<TTree> m_tree;
        std::unique_ptr<TTree> m_run_context_tree;
        ActiveTree m_active_tree = ActiveTree::Event;

        // branch name -> buffer
        std::unordered_map<std::string, std::unique_ptr<IHolder>> m_buffers;
        std::unordered_map<std::string, std::type_index> m_branch_types;

        std::unordered_map<std::string, std::unique_ptr<IHolder>> m_run_context_buffers;
        std::unordered_map<std::string, std::type_index> m_run_context_branch_types;
        std::size_t m_fill_count = 0;
        static constexpr std::size_t m_shrink_interval = 1000;
    };
