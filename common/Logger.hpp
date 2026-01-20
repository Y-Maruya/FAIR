#pragma once

#include <memory>
#include <string>

// spdlog
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/rotating_file_sink.h>

namespace FAIR {

inline std::shared_ptr<spdlog::logger>& default_logger() {
  static std::shared_ptr<spdlog::logger> logger;
  return logger;
}

// Call once at program start (main)
inline void init_logger(const std::string& name = "FAIR",
                        const std::string& logfile = "",
                        spdlog::level::level_enum level = spdlog::level::info) {
  if (default_logger()) return;

  std::shared_ptr<spdlog::logger> lg;

  if (!logfile.empty()) {
    // 20MB x 3 files
    lg = spdlog::rotating_logger_mt(name, logfile, 20 * 1024 * 1024, 3);
  } else {
    lg = spdlog::stdout_color_mt(name);
  }

  lg->set_level(level);
  lg->flush_on(spdlog::level::info);

  // time | level | logger-name | message
  lg->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%^%l%$] [%s:%# %!] %v");

  default_logger() = lg;
  spdlog::set_default_logger(lg);
}

// optional: change log level at runtime
inline void set_level(spdlog::level::level_enum level) {
  if (default_logger()) default_logger()->set_level(level);
  spdlog::set_level(level);
}

} // namespace FAIR

// ---------- Macros (module name is optional) ----------
#define LOG_TRACE(...) SPDLOG_LOGGER_CALL(spdlog::default_logger_raw(), spdlog::level::trace, __VA_ARGS__)
#define LOG_DEBUG(...) SPDLOG_LOGGER_CALL(spdlog::default_logger_raw(), spdlog::level::debug, __VA_ARGS__)
#define LOG_INFO(...)  SPDLOG_LOGGER_CALL(spdlog::default_logger_raw(), spdlog::level::info,  __VA_ARGS__)
#define LOG_WARN(...)  SPDLOG_LOGGER_CALL(spdlog::default_logger_raw(), spdlog::level::warn,  __VA_ARGS__)
#define LOG_ERROR(...) SPDLOG_LOGGER_CALL(spdlog::default_logger_raw(), spdlog::level::err,   __VA_ARGS__)
