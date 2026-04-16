# How to Add a New Reconstruction Algorithm (AHCALRecoAlg)

This document walks through the steps needed to add a new reconstruction or analysis algorithm to the FAIR framework.

## Overview

All algorithms implement the `IAlg` interface (`common/IAlg.hpp`) and are registered in the `AlgRegistry` so they can be instantiated by name from a YAML config file.

The typical layout for a new algorithm `MyAlg` inside `reco_alg/module/` or `analysis/` is:

```
reco_alg/module/MyAlg/
    MyAlg.hpp
    MyAlg.cpp
    CMakeLists.txt
```

---

## Step 1 – Create the header (`MyAlg.hpp`)

```cpp
// reco_alg/module/MyAlg/MyAlg.hpp
#pragma once
#include "common/EventStore.hpp"
#include "common/IAlg.hpp"
#include "common/edm/EDM.hpp"
#include <yaml-cpp/yaml.h>
#include <string>

namespace AHCALRecoAlg {

struct MyAlgCfg {
    std::string in_recohit_key = "RecoHits";   // key to read from EventStore
    std::string out_key        = "MyOutput";   // key to write to EventStore
    double      some_threshold = 10.0;
};

class MyAlg final : public IAlg {
public:
    MyAlg(RunContext& ctx, std::string name)
        : IAlg(ctx, std::move(name)) {}

    void initialize() override;   // optional – called once before processing
    void execute(EventStore& evt) override;
    void finalize() override;     // optional – called once after processing
    void parse_cfg(const YAML::Node& n) override;

private:
    MyAlgCfg m_cfg;
};

} // namespace AHCALRecoAlg
```

**Key points:**
- Place the class in the `AHCALRecoAlg` namespace (conventional for all algs).
- Declare a plain config struct (`MyAlgCfg`) with sensible defaults.
- `execute()` and `parse_cfg()` are pure-virtual in `IAlg` and **must** be implemented.
- `initialize()`, `finalize()`, `init_by_run()`, `finish_by_run()` are optional hooks.

---

## Step 2 – Implement the algorithm (`MyAlg.cpp`)

```cpp
// reco_alg/module/MyAlg/MyAlg.cpp
#include "MyAlg.hpp"
#include "common/AlgRegistry.hpp"
#include "common/config/YAMLUtil.hpp"
#include "common/Logger.hpp"

// Register the algorithm so it can be created by name from YAML.
AHCAL_REGISTER_ALG(AHCALRecoAlg::MyAlg, "MyAlg")

namespace AHCALRecoAlg {

void MyAlg::initialize() {
    LOG_INFO("MyAlg: initialized with threshold = {}", m_cfg.some_threshold);
    // Allocate resources, open files, etc.
}

void MyAlg::execute(EventStore& evt) {
    // Read input from the EventStore.
    auto& hits = evt.get<std::vector<AHCALRecoHit>>(m_cfg.in_recohit_key);
    if (hits.empty()) {
        LOG_DEBUG("MyAlg: no hits in event, skipping.");
        return;
    }

    // Do processing.
    // ...

    // Write output to the EventStore.
    // evt.put(m_cfg.out_key, std::move(result));
}

void MyAlg::finalize() {
    LOG_INFO("MyAlg: finalizing.");
    // Release resources, close files, etc.
}

void MyAlg::parse_cfg(const YAML::Node& n) {
    // get_or<T>(node, "key", default) reads from YAML with a fallback.
    m_cfg.in_recohit_key  = get_or<std::string>(n, "in_recohit_key",  m_cfg.in_recohit_key);
    m_cfg.out_key         = get_or<std::string>(n, "out_key",         m_cfg.out_key);
    m_cfg.some_threshold  = get_or<double>     (n, "some_threshold",  m_cfg.some_threshold);
}

} // namespace AHCALRecoAlg
```

**Key points:**
- Call `AHCAL_REGISTER_ALG(FullyQualifiedType, "StringName")` exactly once per algorithm, at namespace scope in the `.cpp` file.
- Use `evt.get<T>(key)` to read and `evt.put(key, value)` to write data products.
- Use `get_or<T>(node, "key", default)` for YAML config parsing.
- Use `LOG_INFO` / `LOG_DEBUG` / `LOG_WARN` / `LOG_ERROR` for logging.

---

## Step 3 – Create the `CMakeLists.txt`

```cmake
# reco_alg/module/MyAlg/CMakeLists.txt
cmake_minimum_required(VERSION 3.16)

add_library(MyAlg STATIC MyAlg.cpp)

target_include_directories(MyAlg
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${FAIR_ROOT_DIR}
)

if(TARGET fair_options)
  target_link_libraries(MyAlg PUBLIC fair_options)
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(MyAlg PRIVATE -Wall -Wextra -Wpedantic)
endif()

add_library(FAIR::MyAlg ALIAS MyAlg)
```

---

## Step 4 – Register in the parent `CMakeLists.txt`

Open `reco_alg/module/CMakeLists.txt` and add your subdirectory and library:

```cmake
add_subdirectory(MyAlg)   # <-- add this

add_library(reco_module_algs INTERFACE)
target_link_libraries(reco_module_algs INTERFACE
  TrackFitAlg MuonKFAlg RmIsolatedHitAlg SummaryAlg
  MyAlg   # <-- add this
)
```

For an analysis algorithm in `analysis/`, edit `analysis/CMakeLists.txt` the same way.

---

## Step 5 – Add the algorithm to a YAML configuration

`algs` セクションに `MyAlg` を追加します。アルゴリズムは **YAML に記述した順番で実行される**ため、`cfg` で指定する `in_*_key` などの入力キーは、必ず**それより前に並んでいる**アルゴリズムが EventStore に出力済みである必要があります。

```yaml
algs:
  - type: AdcToEnergyReadTTreeAlg   # RecoHits を出力する
    cfg:
      in_rawhit_key: RawHits
      out_recohit_key: RecoHits
      # ... calibration file paths ...

  - type: MyAlg                     # matches the string in AHCAL_REGISTER_ALG
    cfg:
      in_recohit_key: RecoHits      # AdcToEnergyReadTTreeAlg が先に出力していること
      out_key: MyOutput
      some_threshold: 15.0

  - type: RootWriterAlg
    cfg:
      outputlist:
        - SimpleFittedTrack
```

---

## Summary checklist

| Step | What to do |
|------|-----------|
| 1 | Create `reco_alg/module/MyAlg/MyAlg.hpp` – define config struct and `IAlg` subclass |
| 2 | Create `reco_alg/module/MyAlg/MyAlg.cpp` – implement `execute()`, `parse_cfg()`, add `AHCAL_REGISTER_ALG` |
| 3 | Create `reco_alg/module/MyAlg/CMakeLists.txt` – build as a STATIC library |
| 4 | Add `add_subdirectory(MyAlg)` and link in `reco_alg/module/CMakeLists.txt` |
| 5 | Add `MyAlg` to the `algs` list in the YAML config (after any alg that produces its input keys) |
