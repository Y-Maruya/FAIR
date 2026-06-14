#include "common/config/ParseRunConfig.hpp"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <stdexcept>
#include <string>
#include <unistd.h>

namespace {

void expect(bool condition, const std::string& message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

std::string write_snapshot(const std::string& name, const std::string& contents) {
  const std::string path =
      "/tmp/fair_runinfo_snapshot_" + std::to_string(getpid()) + "_" + name + ".json";
  std::ofstream output(path);
  output << contents;
  output.close();
  return path;
}

void test_snapshot_type_and_detail() {
  const std::string path = write_snapshot("valid", R"({
    "version": 1,
    "generated_at": "2026-06-07T00:00:00+00:00",
    "run_types": {
      "22915": "AHCAL",
      "22916": "Physics"
    },
    "details": {
      "22915": {
        "runnumber": 22915,
        "type": "AHCAL",
        "starttime": "2026-05-01T10:00:00.500000 UTC",
        "stoptime": "2026-05-01T10:05:00 UTC",
        "configuration": {
          "components": [{
            "name": "tlureceiver00",
            "modules": [{
              "name": "tlureceiver00",
              "settings": {
                "trigger_logic": "A&B",
                "voltage_threshold": [0.1, 0.2],
                "trigger_stretch": [3, 4],
                "trigger_delay": [1, 2],
                "FASER_CalibRate": 7
              }
            }]
          }]
        }
      }
    }
  })");
  setenv("FAIR_RUNINFO_SNAPSHOT", path.c_str(), 1);

  expect(get_run_type(22915) == "AHCAL", "AHCAL type was not read from snapshot");
  expect(get_run_type(22916) == "Physics", "Physics type was not read from snapshot");
  expect(!find_runinfo_detail_in_snapshot(22916), "Unexpected detail found for run 22916");

  const ConditionStore conditions = parse_condition_store(22915);
  expect(std::abs(conditions.endtime - conditions.starttime - 299.5) < 1e-6,
         "Snapshot timestamps were not parsed correctly");
  expect(conditions.triggerLogic == "A&B", "Trigger logic was not parsed");
  expect(conditions.thresholds == std::vector<double>({0.1, 0.2}),
         "Thresholds were not parsed");
  expect(conditions.triggerStretch == std::vector<int>({3, 4}),
         "Trigger stretch was not parsed");
  expect(conditions.triggerDelay == std::vector<int>({1, 2}),
         "Trigger delay was not parsed");
  expect(conditions.calibRate == 7, "Calibration rate was not parsed");
}

void test_invalid_snapshot_is_unavailable() {
  const std::string path = write_snapshot("invalid", R"({"version": 1, "details": {}})");
  setenv("FAIR_RUNINFO_SNAPSHOT", path.c_str(), 1);

  expect(!find_run_type_in_snapshot(22915).has_value(),
         "Invalid snapshot unexpectedly returned a run type");
  expect(!find_runinfo_detail_in_snapshot(22915),
         "Invalid snapshot unexpectedly returned run details");
}

}  // namespace

int main() {
  test_snapshot_type_and_detail();
  test_invalid_snapshot_is_unavailable();
  return 0;
}
