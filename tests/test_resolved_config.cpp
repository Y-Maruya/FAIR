#include "common/config/ResolvedConfig.hpp"
#include "common/config/YAMLUtil.hpp"

#include <filesystem>
#include <stdexcept>
#include <string>
#include <unistd.h>

namespace {

void expect(bool condition, const std::string& message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

void test_defaults_and_run_values_are_written() {
  YAML::Node config = YAML::Load(R"(
run:
  input: input-*.raw
  output: out.root
  log_file: app.log
  log_level: info
algs:
  - type: ExampleAlg
    cfg:
      supplied: 3
)");

  YAML::Node cfg = config["algs"][0]["cfg"];
  const YAML::Node input_config = YAML::Clone(config);
  expect(get_or<int>(cfg, "supplied", 9) == 3, "Supplied value was not used");
  expect(get_or<double>(cfg, "defaulted", 1.25) == 1.25, "Default value was not used");

  RunConfig run;
  run.input = "/data/run42/*.raw";
  run.output = "/tmp/run42.root";
  run.log_file = "/tmp/fair_resolved_config_test/app.log";
  run.runNumber = 42;
  run.poolIndex = 7;

  const std::filesystem::path output =
      "/tmp/fair_resolved_config_" + std::to_string(getpid()) + ".yaml";
  const std::filesystem::path input_output =
      "/tmp/fair_input_config_" + std::to_string(getpid()) + ".yaml";
  run.config_output = output.string();
  run.config_input_output = input_output.string();
  write_input_config(input_config, run);
  write_resolved_config(config, run);

  const YAML::Node written_input = YAML::LoadFile(input_output.string());
  const YAML::Node written = YAML::LoadFile(output.string());
  expect(!written_input["algs"][0]["cfg"]["defaulted"],
         "Input snapshot unexpectedly contains a parse_cfg default");
  expect(written["run"]["runNumber"].as<int>() == 42, "Run number was not updated");
  expect(written["run"]["input"].as<std::string>() == "/data/run42/*.raw",
         "Run input was not updated");
  expect(written["algs"][0]["cfg"]["supplied"].as<int>() == 3,
         "Supplied algorithm config was not written");
  expect(written["algs"][0]["cfg"]["defaulted"].as<double>() == 1.25,
         "Defaulted algorithm config was not written");

  std::filesystem::remove(input_output);
  std::filesystem::remove(output);
}

void test_default_output_is_per_run() {
  RunConfig first_pool;
  first_pool.log_file = "/tmp/fair/run42/app.log";
  first_pool.runNumber = 42;
  first_pool.poolIndex = 1;

  RunConfig second_pool = first_pool;
  second_pool.poolIndex = 2;

  const std::filesystem::path resolved_path = resolved_config_path(first_pool);
  const std::filesystem::path input_path = input_config_path(first_pool);
  const std::string timestamp = config_execution_timestamp();

  expect(resolved_path ==
             "/tmp/fair/run42/config_run42_" + timestamp + ".yaml",
         "Default resolved config path does not contain the run and execution timestamp");
  expect(resolved_path == resolved_config_path(second_pool),
         "Resolved config path changed between pools of the same run");
  expect(input_path ==
             "/tmp/fair/run42/config_input_run42_" + timestamp + ".yaml",
         "Default input config path does not contain the run and execution timestamp");
  expect(input_path == input_config_path(second_pool),
         "Input config path changed between pools of the same run");
}

}  // namespace

int main() {
  test_defaults_and_run_values_are_written();
  test_default_output_is_per_run();
  return 0;
}
