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
#include <algorithm>

#include "common/Logger.hpp"
#include "common/RunContext.hpp"

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
    load_run_context_entries();
  }

  ~RootInput() { if (m_file) m_file->Close(); }

  Long64_t entries() const { return m_entries; }

  Long64_t current_entry() const { return m_entry - 1; }

  bool read_entry(Long64_t i) {
    if (i < 0 || i >= m_entries) return false;
    m_tree->GetEntry(i);
    m_entry = i + 1;
    return true;
  }

  bool next() {
    if (m_entry >= m_entries) return false;
    LOG_DEBUG("Reading entry {} / {}", m_entry, m_entries);
    m_tree->GetEntry(m_entry++);
    LOG_DEBUG("Read entry {} / {}", m_entry - 1, m_entries);
    return true;
  }

  bool has_branch(const std::string& branch_name) const {
    return m_tree && m_tree->GetBranch(branch_name.c_str());
  }

  bool has_run_context() const { return !m_run_context_entries.empty(); }

  bool apply_run_context_for_event(long long event_counter, RunContext& ctx) const {
    auto it = std::find_if(
      m_run_context_entries.begin(),
      m_run_context_entries.end(),
      [event_counter](const RunContextRange& x) {
        return x.start_event <= event_counter && event_counter <= x.end_event;
      });
    if (it == m_run_context_entries.end()) return false;
    ctx = it->context;
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

  struct RunContextRange {
    long long start_event = 0;
    long long end_event = -1;
    RunContext context;
  };
  std::vector<RunContextRange> m_run_context_entries;

  template <class T>
  static void set_scalar_branch(TTree* tree, const std::string& branch_name, T& value) {
    if (!tree->GetBranch(branch_name.c_str())) {
      throw std::runtime_error("Branch not found: " + branch_name);
    }
    const int rc = tree->SetBranchAddress(branch_name.c_str(), &value);
    if (rc < 0) throw std::runtime_error("SetBranchAddress failed: " + branch_name);
  }

  template <class T>
  static void set_vector_branch(TTree* tree, const std::string& branch_name, std::vector<T>*& value_ptr) {
    if (!tree->GetBranch(branch_name.c_str())) {
      throw std::runtime_error("Branch not found: " + branch_name);
    }
    const int rc = tree->SetBranchAddress(branch_name.c_str(), &value_ptr);
    if (rc < 0) throw std::runtime_error("SetBranchAddress failed: " + branch_name);
  }

  void load_run_context_entries() {
    auto* run_tree = dynamic_cast<TTree*>(m_file->Get("run_context"));
    if (!run_tree) return;

    long long start_event = 0;
    long long end_event = -1;
    std::string input, output, log_file, log_level;
    int run_number = 0;
    int pool_index = 0;
    bool mc = false;
    long long n_events = -1;

    std::vector<int>* skip_layers = nullptr;
    int n_trigger_layers = 0;
    std::vector<int>* trigger_layers = nullptr;
    double starttime = 0.0;
    double endtime = 0.0;
    std::string trigger_logic;
    std::vector<double>* thresholds = nullptr;
    std::vector<int>* trigger_stretch = nullptr;
    std::vector<int>* trigger_delay = nullptr;
    int calib_rate = 0;

    set_scalar_branch(run_tree, "ContextRange.start_event", start_event);
    set_scalar_branch(run_tree, "ContextRange.end_event", end_event);

    set_scalar_branch(run_tree, "RunConfig.input", input);
    set_scalar_branch(run_tree, "RunConfig.output", output);
    set_scalar_branch(run_tree, "RunConfig.log_file", log_file);
    set_scalar_branch(run_tree, "RunConfig.runNumber", run_number);
    set_scalar_branch(run_tree, "RunConfig.poolIndex", pool_index);
    set_scalar_branch(run_tree, "RunConfig.MC", mc);
    set_scalar_branch(run_tree, "RunConfig.nEvents", n_events);
    set_scalar_branch(run_tree, "RunConfig.log_level", log_level);

    set_vector_branch(run_tree, "ConditionStore.skipLayers", skip_layers);
    set_scalar_branch(run_tree, "ConditionStore.nTriggerLayers", n_trigger_layers);
    set_vector_branch(run_tree, "ConditionStore.triggerLayers", trigger_layers);
    set_scalar_branch(run_tree, "ConditionStore.starttime", starttime);
    set_scalar_branch(run_tree, "ConditionStore.endtime", endtime);
    set_scalar_branch(run_tree, "ConditionStore.triggerLogic", trigger_logic);
    set_vector_branch(run_tree, "ConditionStore.thresholds", thresholds);
    set_vector_branch(run_tree, "ConditionStore.triggerStretch", trigger_stretch);
    set_vector_branch(run_tree, "ConditionStore.triggerDelay", trigger_delay);
    set_scalar_branch(run_tree, "ConditionStore.calibRate", calib_rate);

    const Long64_t n_ctx_entries = run_tree->GetEntries();
    m_run_context_entries.reserve(n_ctx_entries);
    for (Long64_t i = 0; i < n_ctx_entries; ++i) {
      run_tree->GetEntry(i);
      RunContextRange r;
      r.start_event = start_event;
      r.end_event = end_event;

      r.context.config.input = input;
      r.context.config.output = output;
      r.context.config.log_file = log_file;
      r.context.config.runNumber = run_number;
      r.context.config.poolIndex = pool_index;
      r.context.config.MC = mc;
      r.context.config.nEvents = n_events;
      r.context.config.log_level = log_level;

      r.context.conditions.skipLayers = skip_layers ? *skip_layers : std::vector<int>{};
      r.context.conditions.nTriggerLayers = n_trigger_layers;
      r.context.conditions.triggerLayers = trigger_layers ? *trigger_layers : std::vector<int>{};
      r.context.conditions.starttime = starttime;
      r.context.conditions.endtime = endtime;
      r.context.conditions.triggerLogic = trigger_logic;
      r.context.conditions.thresholds = thresholds ? *thresholds : std::vector<double>{};
      r.context.conditions.triggerStretch = trigger_stretch ? *trigger_stretch : std::vector<int>{};
      r.context.conditions.triggerDelay = trigger_delay ? *trigger_delay : std::vector<int>{};
      r.context.conditions.calibRate = calib_rate;

      m_run_context_entries.push_back(std::move(r));
    }
    LOG_INFO("Loaded {} run_context entries from input ROOT file", m_run_context_entries.size());
  }
};
