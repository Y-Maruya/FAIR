// IO/Descriptor.hpp
#pragma once
#include "IO/writer/RootOutput.hpp"
#include "IO/reader/RootInput.hpp"
#include <string>
#include <vector>
#include <functional>
#include <type_traits>

struct FieldDesc {
  std::string name;
  std::function<void(const void* obj, RootOutput& out, const std::string& prefix)> write;
  std::function<void(void* obj, RootInput& in, const std::string& prefix)> read;
};

template <class T, class M>
FieldDesc field(std::string n, M T::*member) {
  FieldDesc d;
  d.name = std::move(n);
  const std::string name_copy = d.name;

  d.write = [member, name_copy](const void* obj, RootOutput& out, const std::string& prefix) {
    const T& x = *static_cast<const T*>(obj);
    using MT = std::decay_t<decltype(x.*member)>;
    *out.get_or_make_buffer<MT>(prefix + "." + name_copy) = x.*member;
  };

  d.read = [member, name_copy](void* obj, RootInput& in, const std::string& prefix) {
    T& x = *static_cast<T*>(obj);
    using MT = std::decay_t<decltype(x.*member)>;
    const MT* buf = in.get_or_make_address<MT>(prefix + "." + name_copy);
    x.*member = *buf;
  };

  return d;
}

struct FieldDescVector {
  std::string name;
  std::function<void(const std::vector<const void*>& objs, RootOutput& out, const std::string& prefix)> write;

  std::function<std::size_t(RootInput& in, const std::string& prefix)> size;
  std::function<void(const std::vector<void*>& objs, RootInput& in, const std::string& prefix)> read;
};

template <class T, class M>
FieldDescVector field_vector(std::string n, M T::*member) {
  FieldDescVector d;
  d.name = std::move(n);
  const std::string name_copy = d.name;

  d.write =
    [member, name_copy](const std::vector<const void*>& objs, RootOutput& out, const std::string& prefix) {
      using MT = std::decay_t<M>;
      auto* buf = out.get_or_make_buffer<std::vector<MT>>(prefix + "." + name_copy);
      buf->clear();
      buf->reserve(objs.size());
      for (const void* p : objs) {
        const T& x = *static_cast<const T*>(p);
        buf->push_back(x.*member);
      }
    };

  d.size = [name_copy](RootInput& in, const std::string& prefix) -> std::size_t {
    using MT = std::decay_t<M>;
    const auto* buf = in.get_or_make_address<std::vector<MT>>(prefix + "." + name_copy);
    LOG_DEBUG("Getting size of vector field '{}': {}", prefix + "." + name_copy, buf->size());
    return buf->size();
  };

  d.read =
    [member, name_copy](const std::vector<void*>& objs, RootInput& in, const std::string& prefix) {
      using MT = std::decay_t<M>;
      const auto* buf = in.get_or_make_address<std::vector<MT>>(prefix + "." + name_copy);
      const std::size_t n = std::min(objs.size(), buf->size());
      for (std::size_t i = 0; i < n; ++i) {
        T& x = *static_cast<T*>(objs[i]);
        x.*member = (*buf)[i];
      }
    };

  return d;
}