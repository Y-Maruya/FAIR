#pragma once
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>

#include <memory>
#include <string>
#include <typeindex>
#include <unordered_map>
#include <stdexcept>
#include <vector>
#include <type_traits>

#include "common/Logger.hpp"

class RootInput {
public:
  RootInput(const std::string& filename, const std::string& treename = "events")
    : m_file(TFile::Open(filename.c_str(), "READ")) {
    if (!m_file || m_file->IsZombie()) {
      LOG_ERROR("Failed to open ROOT file: {}", filename);
      throw std::runtime_error("Failed to open ROOT file: " + filename);
    }
    m_tree = dynamic_cast<TTree*>(m_file->Get(treename.c_str()));
    if (!m_tree) throw std::runtime_error("TTree not found: " + treename);
    m_entries = m_tree->GetEntries();
  }

  ~RootInput() { if (m_file) m_file->Close(); }

  Long64_t entries() const { return m_entries; }

  // 直近に GetEntry したエントリ番号（未読なら -1）
  Long64_t current_entry() const { return m_entry - 1; }

  bool read_entry(Long64_t i) {
    if (i < 0 || i >= m_entries) return false;
    m_tree->GetEntry(i);
    m_entry = i + 1;
    return true;
  }

  bool next() {
    if (m_entry >= m_entries) return false;
    m_tree->GetEntry(m_entry++);
    return true;
  }

  template <class T>
  T* get_or_make_address(const std::string& branch_name) {
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
      holder->set_branch(m_tree, branch_name);

      T* raw = holder->ptr();
      m_branch_types.emplace(branch_name, want);
      m_buffers.emplace(branch_name, std::move(holder));
      if (m_entry > 0) {
        const Long64_t nb = m_tree->GetEntry(m_entry - 1);
        LOG_DEBUG("Reload entry {} after binding '{}': bytes={}", m_entry - 1, branch_name, nb);
      }
      return raw;
    }

    return static_cast<Holder<T>*>(it->second.get())->ptr();
  }

private:
  template <typename>
  struct is_std_vector : std::false_type {};
  template <typename U, typename A>
  struct is_std_vector<std::vector<U, A>> : std::true_type {};

  struct IHolder { virtual ~IHolder() = default; };

  template <class T>
  struct Holder final : IHolder {
    T value{};          // non-vector
    T* vec_ptr = nullptr; // vector only

    void set_branch(TTree* t, const std::string& bn) {
      if (!t) throw std::runtime_error("TTree is null");
      if (!t->GetBranch(bn.c_str())) {
        LOG_ERROR("Branch not found: '{}'", bn);
        LOG_ERROR("Please check the name");
        throw std::runtime_error("Branch not found: " + bn);
      }

      int rc = 0;
      if constexpr (is_std_vector<T>::value) {
        if (!vec_ptr) vec_ptr = new T();
        rc = t->SetBranchAddress(bn.c_str(), &vec_ptr);
      } else {
        rc = t->SetBranchAddress(bn.c_str(), &value);
      }

      if (rc < 0) {
        LOG_ERROR("SetBranchAddress failed for '{}' (rc={}) type={}", bn, rc, typeid(T).name());
        throw std::runtime_error("SetBranchAddress failed: " + bn);
      }

      LOG_DEBUG("SetBranchAddress '{}' type={}", bn, typeid(T).name());
    }

    ~Holder() override {
      if constexpr (is_std_vector<T>::value) {
        delete vec_ptr;
        vec_ptr = nullptr;
      }
    }

    T* ptr() {
      if constexpr (is_std_vector<T>::value) return vec_ptr;
      return &value;
    }
  };

  std::unique_ptr<TFile> m_file;
  TTree* m_tree = nullptr;
  Long64_t m_entries = 0;
  Long64_t m_entry = 0;

  std::unordered_map<std::string, std::unique_ptr<IHolder>> m_buffers;
  std::unordered_map<std::string, std::type_index> m_branch_types;
};