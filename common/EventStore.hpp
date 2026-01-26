// event_store.hpp
#pragma once
#include <any>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <string_view>
#include <typeindex>
#include <unordered_map>
#include <utility>
#include <memory>
#include <vector>
#include "common/Logger.hpp"
// A tiny, type-safe-ish event store: per-event key-value container.
// - EventStore definition never changes when you add new RecoAlg outputs.
// - Stores objects by value (inside std::any). Use move to avoid copies.
// - Type mismatch / missing key throws with helpful error.

class EventStore {
public:
  // Put (copy)
  template <class T>
  void put(std::string key, const T& value) {
    Item item;
    item.type = std::type_index(typeid(T));
    item.payload = value;  // copy
    m_map.emplace(std::move(key), std::move(item));
  }

  // Put (move)
  template <class T>
  void put(std::string key, T&& value) {
    Item item;
    item.type = std::type_index(typeid(T));
    item.payload = std::forward<T>(value);  // move if possible
    m_map.emplace(std::move(key), std::move(item));
  }

  // Overwrite existing (copy/move)
  template <class T>
  void set(std::string key, T&& value) {
    Item item;
    item.type = std::type_index(typeid(T));
    item.payload = std::forward<T>(value);
    m_map[std::move(key)] = std::move(item);
  }

  bool has(std::string_view key) const {
    return m_map.find(std::string(key)) != m_map.end();
  }

  void erase(std::string_view key) {
    m_map.erase(std::string(key));
  }
  // Return list of keys (copy). Cheap enough; typical keys count is small.
  std::vector<std::string> keys() const {
    std::vector<std::string> ks;
    ks.reserve(m_map.size());
    for (const auto& kv : m_map){
      ks.push_back(kv.first);
      LOG_DEBUG("EventStore key: '{}'", kv.first);
    }
    return ks;
  }

  // Access stored payload as std::any (const)
  const std::any& any(std::string_view key) const {
    auto it = m_map.find(std::string(key));
    if (it == m_map.end()) {
      LOG_ERROR("EventStore::any: missing key '{}'", key);
      throw std::runtime_error("EventStore::any: missing key '" + std::string(key) + "'");
    }
    return it->second.payload;
  }

  // Optional: non-const any() (rarely needed, but symmetrical)
  std::any& any(std::string_view key) {
    auto it = m_map.find(std::string(key));
    if (it == m_map.end()) {
      LOG_ERROR("EventStore::any: missing key '{}'", key);
      throw std::runtime_error("EventStore::any: missing key '" + std::string(key) + "'");
    }
    return it->second.payload;
  }
  void clear() { m_map.clear(); }

  // Get mutable reference
  template <class T>
  T& get(std::string_view key) {
    auto it = m_map.find(std::string(key));
    if (it == m_map.end()) {
      LOG_ERROR("EventStore::get: missing key '{}'", key);
      throw std::runtime_error("EventStore::get: missing key '" + std::string(key) + "'");
    }
    if (it->second.type != std::type_index(typeid(T))) {
      LOG_ERROR(
        "EventStore::get: type mismatch for key '{}' (stored={}, requested={})",
        key, it->second.type.name(), typeid(T).name()
      );
      throw std::runtime_error(
        "EventStore::get: type mismatch for key '" + std::string(key) +
        "' (stored=" + std::string(it->second.type.name()) +
        ", requested=" + std::string(typeid(T).name()) + ")"
      );
    }
    return std::any_cast<T&>(it->second.payload);
  }

  // Get const reference
  template <class T>
  const T& get(std::string_view key) const {
    auto it = m_map.find(std::string(key));
    if (it == m_map.end()) {
      LOG_ERROR("EventStore::get: missing key '{}'", key);
      throw std::runtime_error("EventStore::get: missing key '" + std::string(key) + "'");
    }
    if (it->second.type != std::type_index(typeid(T))) {
      LOG_ERROR(
        "EventStore::get: type mismatch for key '{}' (stored={}, requested={})",
        key, it->second.type.name(), typeid(T).name()
      );
      throw std::runtime_error(
        "EventStore::get: type mismatch for key '" + std::string(key) +
        "' (stored=" + std::string(it->second.type.name()) +
        ", requested=" + std::string(typeid(T).name()) + ")"
      );
    }
    return std::any_cast<const T&>(it->second.payload);
  }

  // Optional get (returns nullptr if missing or type mismatch)
  template <class T>
  T* try_get(std::string_view key) {
    auto it = m_map.find(std::string(key));
    if (it == m_map.end()) return nullptr;
    if (it->second.type != std::type_index(typeid(T))) return nullptr;
    return &std::any_cast<T&>(it->second.payload);
  }

private:
  struct Item {
    std::type_index type{typeid(void)};
    std::any payload;
  };

  std::unordered_map<std::string, Item> m_map;
};
