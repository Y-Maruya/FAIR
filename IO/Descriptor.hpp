// IO/Descriptor.hpp
#pragma once
#include "IO/writer/RootOutput.hpp"
#include <string>
#include <vector>
#include <functional>
#include <type_traits>

// class RootOutput; // forward

struct FieldDesc {
  std::string name;
  std::function<void(const void* obj, RootOutput& out, const std::string& prefix)> write;
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
  return d;
}

struct FieldDescVector {
    std::string name;
    std::function<void(const std::vector<const void*>& objs, RootOutput& out, const std::string& prefix)> write;
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
  return d;
}