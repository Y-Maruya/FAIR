# How to Save Data to ROOT

This document explains two approaches to writing analysis results to ROOT files in the FAIR framework:

1. **Framework-managed output** – register your EDM struct and let `RootWriterAlg` persist it automatically.
2. **Self-managed ROOT output** – open and write your own `TFile` / `TTree` / histograms inside an analysis algorithm (e.g. `MuonEffAlg`).

---

## Approach 1 – Framework-managed output via `RootWriterAlg`

Use this when you want to persist a structured data product (a plain C++ struct or a vector of structs) to a ROOT TTree so it can be read back later by another algorithm or processed in a ROOT macro.

### 1-a Define your EDM struct

Add a new header under `common/edm/` (or alongside your algorithm):

```cpp
// common/edm/MyResult.hpp
#pragma once
#include "IO/Descriptor.hpp"
#include "IO/IOTypeRegistry.hpp"

struct MyResult {
    double value = 0.0;
    int    count = 0;
    bool   valid = false;
};

// Describe the fields that should be written/read.
// The first argument to field() is the branch name inside the TTree.
inline std::vector<FieldDesc> describe(const MyResult*) {
    return {
        field("value", &MyResult::value),
        field("count", &MyResult::count),
        field("valid", &MyResult::valid),
    };
}

// Register for both writer and reader so the type is known to the framework.
AHCAL_REGISTER_IO_STRUCT(MyResult, "MyResult");
```

For a `vector<MyResult>` you also need:

```cpp
inline std::vector<FieldDescVector> describe_vector(const MyResult*) {
    return {
        field_vector("v.value", &MyResult::value),
        field_vector("v.count", &MyResult::count),
        field_vector("v.valid", &MyResult::valid),
    };
}

AHCAL_REGISTER_IO_STRUCT_VECTOR(MyResult, "vector<MyResult>");
```

**Reference:** `common/edm/RecoHit.hpp` shows both macros in use.

### 1-b Include the header in `common/edm/EDM.hpp`

```cpp
// common/edm/EDM.hpp
#pragma once
#include "common/edm/RecoHit.hpp"
// ...existing includes...
#include "common/edm/MyResult.hpp"   // <-- add this
```

### 1-c Put the result in the EventStore from your algorithm

```cpp
void MyAlg::execute(EventStore& evt) {
    MyResult result;
    result.value = 42.0;
    result.count = 7;
    result.valid = true;
    evt.put(m_cfg.out_key, std::move(result));
}
```

### 1-d Declare it in the YAML `outputlist`

```yaml
algs:
  # ... your processing algs ...
  - type: RootWriterAlg
    cfg:
      outputlist:
        - MyResult          # single struct
        # - vector<MyResult>  # if you used AHCAL_REGISTER_IO_STRUCT_VECTOR
```

The framework will create a branch per field under a TTree named after the type string. The output file is the path set in `run.output`.

### 1-e Reading back with `RootInput`

前の処理で書き出したROOTファイルを後続の処理で読み込む場合、`run.input` にROOTファイルのパス（またはglobパターン）を、`reader.type` に `RootInput` を指定します。

```yaml
run:
  input: /path/to/output/my_alg_out.root
  output: /path/to/downstream_out.root
  log_file: /path/to/log
  nEvents: -1
  log_level: INFO
  MC: false
  runNumber: 21723
  poolIndex: 0

reader:
  type: RootInput
  cfg:
    inputlist:
      - [MyResult, MyResultKey]
```

---

## Approach 2 – Self-managed ROOT output inside an analysis algorithm

Use this when you need full control over histograms, `TEfficiency`, `TCanvas`, custom `TTree` layout, etc.

**Reference implementation:** `analysis/MuonEffAlg/MuonEffAlg.cpp` / `MuonEffAlg.hpp`.

### 2-a Store a `TFile` and histogram pointers in `Impl`

The `Impl` pimpl idiom keeps ROOT headers out of the public `.hpp`:

```cpp
// MyAnaAlg.hpp
#pragma once
#include "common/IAlg.hpp"
#include "common/EventStore.hpp"
#include <yaml-cpp/yaml.h>
#include <memory>
#include <string>

namespace AHCALRecoAlg {

struct MyAnaAlgCfg {
    std::string in_track_key  = "FittedTrack";
    std::string out_filename  = "my_ana_out.root";
};

class MyAnaAlg final : public IAlg {
public:
    MyAnaAlg(RunContext& ctx, std::string name)
        : IAlg(ctx, std::move(name)) {}
    ~MyAnaAlg() override;

    void initialize() override;
    void execute(EventStore& evt) override;
    void finalize() override;
    void parse_cfg(const YAML::Node& cfg) override;

private:
    MyAnaAlgCfg cfg_;

    // Pimpl – keeps TFile/TH1 types out of the header.
    struct Impl;
    struct ImplDeleter { void operator()(Impl* p) const; };
    std::unique_ptr<Impl, ImplDeleter> impl_;
};

} // namespace AHCALRecoAlg
```

### 2-b Implement `Impl` and the lifecycle methods in `.cpp`

```cpp
// MyAnaAlg.cpp
#include "MyAnaAlg.hpp"
#include "common/AlgRegistry.hpp"
#include "common/config/YAMLUtil.hpp"
#include "common/Logger.hpp"
#include "common/edm/SimpleFittedTrack.hpp"

#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>

AHCAL_REGISTER_ALG(AHCALRecoAlg::MyAnaAlg, "MyAnaAlg")

namespace AHCALRecoAlg {

// ---- Impl definition (ROOT objects live here) ----
struct MyAnaAlg::Impl {
    explicit Impl(const MyAnaAlgCfg& cfg) {
        fout.reset(TFile::Open(cfg.out_filename.c_str(), "RECREATE"));
        if (!fout || fout->IsZombie())
            throw std::runtime_error("MyAnaAlg: cannot open output file: " + cfg.out_filename);

        h_chi2 = std::make_unique<TH1D>("h_chi2", "Track #chi^{2}/ndf;#chi^{2}/ndf;entries", 100, 0, 20);
    }

    void fill(const SimpleFittedTrack& t) {
        if (!t.valid) return;
        h_chi2->Fill(t.chi2xperndf);
    }

    void write() {
        fout->cd();
        h_chi2->Write();
        LOG_INFO("MyAnaAlg: wrote output to {}", fout->GetName());
    }

    std::unique_ptr<TFile> fout;
    std::unique_ptr<TH1D>  h_chi2;
};

void MyAnaAlg::ImplDeleter::operator()(Impl* p) const { delete p; }

// ---- IAlg lifecycle ----
void MyAnaAlg::initialize() {
    impl_.reset(new Impl(cfg_));
    LOG_INFO("MyAnaAlg: initialized");
}

void MyAnaAlg::execute(EventStore& evt) {
    auto* track = evt.try_get<SimpleFittedTrack>(cfg_.in_track_key);
    if (!track) {
        LOG_DEBUG("MyAnaAlg: no track in event.");
        return;
    }
    impl_->fill(*track);
}

void MyAnaAlg::finalize() {
    if (impl_) impl_->write();
}

MyAnaAlg::~MyAnaAlg() = default;

void MyAnaAlg::parse_cfg(const YAML::Node& n) {
    cfg_.in_track_key = get_or<std::string>(n, "in_track_key", cfg_.in_track_key);
    cfg_.out_filename = get_or<std::string>(n, "out_filename", cfg_.out_filename);
}

} // namespace AHCALRecoAlg
```

**Key points:**
- Open `TFile` in `initialize()`, not in the constructor, so the `RunContext` is fully set up.
- Fill histograms/trees inside `execute()` per event.
- Call `Write()` and close the file in `finalize()`, which runs once after all events are processed.
- Use `evt.try_get<T>(key)` instead of `evt.get<T>(key)` when the product may be absent.

### 2-c `CMakeLists.txt` for the analysis algorithm

```cmake
# analysis/MyAnaAlg/CMakeLists.txt
cmake_minimum_required(VERSION 3.16)

add_library(MyAnaAlg STATIC MyAnaAlg.cpp)

target_include_directories(MyAnaAlg PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${FAIR_ROOT_DIR}
)

if(TARGET fair_options)
  target_link_libraries(MyAnaAlg PUBLIC fair_options)
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(MyAnaAlg PRIVATE -Wall -Wextra -Wpedantic)
endif()

add_library(FAIR::MyAnaAlg ALIAS MyAnaAlg)
```

### 2-d Register in `analysis/CMakeLists.txt`

```cmake
add_subdirectory(MyAnaAlg)   # <-- add

add_library(analysis_algs INTERFACE)
target_link_libraries(analysis_algs INTERFACE
  MuonEffAlg VetoEffAlg InputEffAlg NoiseHitAlg RateAnaAlg VetoAnaAlg
  MyAnaAlg   # <-- add
)
```

### 2-e Include in `AlgFactory.hpp`

```cpp
// common/AlgFactory.hpp
#include "analysis/MyAnaAlg/MyAnaAlg.hpp"   // <-- add
```

### 2-f YAML configuration

`run.input` には生データ（`.raw`）へのパスまたはglobパターンを指定します。

```yaml
run:
  input: /eos/experiment/faser/raw/2026/021723/FaserAHCAL-Physics-021723-*.raw
  output: /path/to/output/my_run_out.root
  log_file: /path/to/output/log
  nEvents: -1
  log_level: INFO
  MC: false
  runNumber: 21723
  poolIndex: 0

reader:
  type: BinaryRawHitReader
  cfg:
    out_rawhits_key: RawHits
    out_tlu_key: TLURawData

algs:
  - type: AdcToEnergyReadTTreeAlg
    cfg:
      in_rawhit_key: RawHits
      out_recohit_key: RecoHits
      # ... calibration file paths ...

  - type: TrackFitAlg
    cfg:
      in_recohit_key: RecoHits
      out_track_key: FittedTrack

  - type: MyAnaAlg
    cfg:
      in_track_key: FittedTrack
      out_filename: /path/to/output/my_ana_out.root
```

---

## Comparison

| | Framework-managed (`RootWriterAlg`) | Self-managed (`TFile` in algorithm) |
|---|---|---|
| Setup effort | Low – describe fields + register | Higher – manage TFile / histograms manually |
| Output format | TTree with one branch per field | Fully custom (histograms, trees, canvases, …) |
| Read back in FAIR | Yes – via `RootInput` | Only if you implement your own reader |
| Best for | Storing structured data products for downstream processing | Producing final analysis plots and summaries |
