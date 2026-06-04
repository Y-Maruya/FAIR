#include "common/config/ExpandConfig.hpp"

#include <cstdlib>
#include <stdexcept>
#include <string>

namespace {

class ScopedEnv {
 public:
  ScopedEnv(const char* name, const char* value) : name_(name) {
    const char* old_value = std::getenv(name_.c_str());
    if (old_value != nullptr) {
      had_old_value_ = true;
      old_value_ = old_value;
    }
    setenv(name_.c_str(), value, 1);
  }

  ~ScopedEnv() {
    if (had_old_value_) {
      setenv(name_.c_str(), old_value_.c_str(), 1);
    } else {
      unsetenv(name_.c_str());
    }
  }

 private:
  std::string name_;
  bool had_old_value_ = false;
  std::string old_value_;
};

void expect(bool condition, const std::string& message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

void test_expand_string() {
  ScopedEnv env("FAIR_BASE", "/tmp/fair");
  expect(expand_env_string("${FAIR_BASE}/calibration/reffiles/mip.root") ==
             "/tmp/fair/calibration/reffiles/mip.root",
         "FAIR_BASE was not expanded in a string");
}

void test_missing_env_throws() {
  unsetenv("FAIR_TEST_UNSET_DO_NOT_DEFINE");

  bool threw = false;
  try {
    (void)expand_env_string("${FAIR_TEST_UNSET_DO_NOT_DEFINE}/x.root");
  } catch (const std::runtime_error& ex) {
    threw = std::string(ex.what()).find("FAIR_TEST_UNSET_DO_NOT_DEFINE") != std::string::npos;
  }

  expect(threw, "Missing environment variable did not throw");
}

void test_yaml_traversal_and_scalar_types() {
  ScopedEnv env("FAIR_BASE", "/tmp/fair");
  YAML::Node node = YAML::Load(R"(
run:
  input: ${FAIR_BASE}/input.raw
  MC: false
  nEvents: -1
algs:
  - cfg:
      file: ${FAIR_BASE}/calibration/reffiles/dac.root
)");

  expand_env_in_yaml(node);

  expect(node["run"]["input"].as<std::string>() == "/tmp/fair/input.raw",
         "Nested run input was not expanded");
  expect(node["algs"][0]["cfg"]["file"].as<std::string>() ==
             "/tmp/fair/calibration/reffiles/dac.root",
         "Nested sequence value was not expanded");
  expect(node["run"]["MC"].as<bool>() == false, "Boolean scalar changed unexpectedly");
  expect(node["run"]["nEvents"].as<int>() == -1, "Integer scalar changed unexpectedly");
}

}  // namespace

int main() {
  test_expand_string();
  test_missing_env_throws();
  test_yaml_traversal_and_scalar_types();
  return 0;
}
