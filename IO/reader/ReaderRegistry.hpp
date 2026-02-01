#pragma once
#include <any>
#include <functional>
#include <string>
#include <unordered_map>
#include <vector>

#include "IO/Descriptor.hpp"
#include "IO/reader/RootInput.hpp"

class ReaderRegistry {
public:
  using ReaderFn = std::function<std::any(const std::string& prefix, RootInput& in)>;

  template <class T>
  void register_struct(std::string type_name) {
    m_readers.emplace(std::move(type_name),
      [](const std::string& prefix, RootInput& in) -> std::any {
        T obj{};
        const auto& desc = describe((const T*)nullptr);
        for (const auto& f : desc) f.read(&obj, in, prefix);
        return obj;
      });
  }

  template <class T>
  void register_vector_struct(std::string type_name) {
    m_readers.emplace(std::move(type_name),
      [](const std::string& prefix, RootInput& in) -> std::any {
        std::vector<T> vec;
        const auto& desc = describe_vector((const T*)nullptr);
        if (desc.empty()) return vec;
        const std::size_t n = desc.front().size(in, prefix);
        LOG_DEBUG("Vector size of desc.front '{}': {}", desc.front().name, n);
        vec.resize(n);

        std::vector<void*> ptrs;
        ptrs.reserve(n);
        for (auto& x : vec) ptrs.push_back(&x);

        for (const auto& f : desc) f.read(ptrs, in, prefix);
        return vec;
      });
  }

  std::any read_any(const std::string& type_name, const std::string& prefix, RootInput& in) const {
    return m_readers.at(type_name)(prefix, in);
  }
  template <class T>
  T read(const std::string& type_name, const std::string& prefix, RootInput& in) const {
    auto any = read_any(type_name, prefix, in);
    return std::any_cast<T>(any);
  }

private:
  std::unordered_map<std::string, ReaderFn> m_readers;
};