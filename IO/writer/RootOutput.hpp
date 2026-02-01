#pragma once
#include <TFile.h>
#include <TTree.h>
#include "common/Logger.hpp"
#include <memory>
#include <string>
#include <unordered_map>
#include <typeindex>
#include <stdexcept>


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
    }

    ~RootOutput() {
        if (m_file) {
            m_file->cd();
            m_tree->Write();

            // Detach before closing so the file doesn't try to own/delete it.
            m_tree->SetDirectory(nullptr);

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
        const std::type_index want(typeid(T));

        auto it_type = m_branch_types.find(branch_name);
        if (it_type != m_branch_types.end() && it_type->second != want) {
        LOG_ERROR("Branch '{}' requested with different type. existing={}, requested={}",
                    branch_name, it_type->second.name(), want.name());
        throw std::runtime_error("Branch type mismatch: " + branch_name);
        }

        auto it = m_buffers.find(branch_name);
        if (it == m_buffers.end()) {
            auto holder = std::make_unique<Holder<T>>();
            T* raw_ptr = &holder->value;

            m_tree->Branch(branch_name.c_str(), raw_ptr);
            LOG_DEBUG("Created branch '{}' of type {}", branch_name, want.name());
            m_branch_types.emplace(branch_name, want);
            m_buffers.emplace(branch_name, std::move(holder));
            return raw_ptr;
        }

        return static_cast<Holder<T>*>(it->second.get())->ptr();
    }

    void fill() {
        m_tree->Fill();
    }

    private:
        struct IHolder {
            virtual ~IHolder() = default;
        };

        template <class T>
        struct Holder final : IHolder {
            T value{};
            T* ptr() { return &value; }
        };

        std::unique_ptr<TFile> m_file;
        std::unique_ptr<TTree> m_tree;

        // branch name -> buffer
        std::unordered_map<std::string, std::unique_ptr<IHolder>> m_buffers;
        std::unordered_map<std::string, std::type_index> m_branch_types;
    };
