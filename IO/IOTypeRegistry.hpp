#pragma once
#include <string>
#include <unordered_map>

#include "common/EventStore.hpp"
#include "IO/reader/ReaderRegistry.hpp"
#include "IO/writer/WriterRegistry.hpp"

// ---- macro utilities (unique internal symbol names) ----
#define AHCAL_DETAIL_CONCAT_INNER(a, b) a##b
#define AHCAL_DETAIL_CONCAT(a, b) AHCAL_DETAIL_CONCAT_INNER(a, b)
#define AHCAL_DETAIL_UNIQUE_NAME(prefix) AHCAL_DETAIL_CONCAT(prefix, __COUNTER__)

class IOTypeRegistry {
public:
  using WriterFn  = void (*)(WriterRegistry&);
  using ReaderFn  = void (*)(ReaderRegistry&, const std::string& type_name);
  using ReadPutFn = void (*)(EventStore&, ReaderRegistry&, RootInput&, const std::string& type_name, const std::string& key);

  static IOTypeRegistry& instance() {
    static IOTypeRegistry inst;
    return inst;
  }

  void add_writer(const std::string& name, WriterFn fn)   { writer_[name]  = fn; }
  void add_reader(const std::string& name, ReaderFn fn)   { reader_[name]  = fn; }
  void add_readput(const std::string& name, ReadPutFn fn) { readput_[name] = fn; }

  WriterFn get_writer(const std::string& name) const {
    auto it = writer_.find(name);
    return (it == writer_.end()) ? nullptr : it->second;
  }
  ReaderFn get_reader(const std::string& name) const {
    auto it = reader_.find(name);
    return (it == reader_.end()) ? nullptr : it->second;
  }
  ReadPutFn get_readput(const std::string& name) const {
    auto it = readput_.find(name);
    return (it == readput_.end()) ? nullptr : it->second;
  }

private:
  std::unordered_map<std::string, WriterFn>  writer_;
  std::unordered_map<std::string, ReaderFn>  reader_;
  std::unordered_map<std::string, ReadPutFn> readput_;
};

namespace io_registry_detail {

// writer
template <class T>
inline void register_struct_writer(WriterRegistry& reg) {
  reg.register_struct<T>();
}
template <class ElemT>
inline void register_vector_elem_writer(WriterRegistry& reg) {
  reg.register_vector_struct<ElemT>();
}

// reader (NOTE: ReaderRegistry wants type_name, not key)
template <class T>
inline void register_struct_reader(ReaderRegistry& reg, const std::string& type_name) {
  reg.register_struct<T>(type_name);
}
template <class ElemT>
inline void register_vector_elem_reader(ReaderRegistry& reg, const std::string& type_name) {
  reg.register_vector_struct<ElemT>(type_name);
}

// read + put
template <class T>
inline void readput_struct(EventStore& store, ReaderRegistry& rr, RootInput& in,
                           const std::string& type_name, const std::string& key) {
  auto obj = rr.read<T>(type_name, key, in);
  store.put(key, std::move(obj));
}
template <class ElemT>
inline void readput_vector(EventStore& store, ReaderRegistry& rr, RootInput& in,
                           const std::string& type_name, const std::string& key) {
  auto vec = rr.read<std::vector<ElemT>>(type_name, key, in);
  store.put(key, std::move(vec));
}

} // namespace io_registry_detail

// Header-only self-registration helpers (place in EDM headers)
#define AHCAL_REGISTER_IO_STRUCT(Type, NameStr)                                      \
  namespace {                                                                        \
  inline const bool AHCAL_DETAIL_UNIQUE_NAME(ahcal_io_reg_) = []() {                 \
    IOTypeRegistry::instance().add_writer(NameStr, &io_registry_detail::register_struct_writer<Type>); \
    IOTypeRegistry::instance().add_reader(NameStr, &io_registry_detail::register_struct_reader<Type>); \
    IOTypeRegistry::instance().add_readput(NameStr, &io_registry_detail::readput_struct<Type>);        \
    return true;                                                                     \
  }();                                                                               \
  }

#define AHCAL_REGISTER_IO_VECTOR_ELEM(ElemType, NameStr)                             \
  namespace {                                                                        \
  inline const bool AHCAL_DETAIL_UNIQUE_NAME(ahcal_io_reg_) = []() {                 \
    IOTypeRegistry::instance().add_writer(NameStr, &io_registry_detail::register_vector_elem_writer<ElemType>); \
    IOTypeRegistry::instance().add_reader(NameStr, &io_registry_detail::register_vector_elem_reader<ElemType>); \
    IOTypeRegistry::instance().add_readput(NameStr, &io_registry_detail::readput_vector<ElemType>);              \
    return true;                                                                     \
  }();                                                                               \
  }

// Backward-compat alias (your EDM headers already use this name)
#define AHCAL_REGISTER_IO_STRUCT_VECTOR(ElemType, NameStr) AHCAL_REGISTER_IO_VECTOR_ELEM(ElemType, NameStr)