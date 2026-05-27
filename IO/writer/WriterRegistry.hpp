// IO/WriterRegistry.hpp
#pragma once
#include <any>
#include <typeindex>
#include <unordered_map>
#include <unordered_set>
#include <functional>
#include <string>
#include "IO/writer/RootOutput.hpp"
#include "IO/Descriptor.hpp"
#include "common/Logger.hpp"

class WriterRegistry {
public:
    using WriterFn = std::function<void(const std::string& key, const std::any& a, RootOutput& out)>;

    template <class T>
    void register_struct() {
        std::type_index ti(typeid(T));
        m_writers[ti] = [](const std::string& key, const std::any& a, RootOutput& out) {
            const T& obj = std::any_cast<const T&>(a);
            const auto& desc = describe((const T*)nullptr);
            for (const auto& f : desc) f.write(&obj, out, key);
        };
        LOG_DEBUG("Registered writer for type '{}'", ti.name());
    }

    template<class T>
    void register_vector_struct() {
        using X = std::vector<T>;
        std::type_index ti(typeid(X));
        if (m_writers.count(ti)) return;

        m_writers[ti] = [](const std::string& key, const std::any& a, RootOutput& out) {
            const auto& vec = std::any_cast<const std::vector<T>&>(a);
            const auto& desc = describe_vector((const T*)nullptr);

            std::vector<const void*> ptrs;
            ptrs.reserve(vec.size());
            for (const auto& x : vec) ptrs.push_back(&x);

            for (const auto& f : desc) f.write(ptrs, out, key);
        };
    }

    bool can_write(const std::any& a) const {
        return m_writers.find(std::type_index(a.type())) != m_writers.end();
    }

    bool is_key_enabled(const std::string& key) const {
        // If whitelist is empty, all keys are enabled (backward compatible)
        if (m_key_whitelist.empty()) {
            return true;
        }
        return m_key_whitelist.find(key) != m_key_whitelist.end();
    }

    void write_any(const std::string& key, const std::any& a, RootOutput& out) const {
        if (is_key_enabled(key)) {
            m_writers.at(std::type_index(a.type()))(key, a, out);
        }
    }

    void set_key_whitelist(const std::unordered_set<std::string>& keys) {
        m_key_whitelist = keys;
        LOG_DEBUG("Set WriterRegistry key whitelist with {} keys", keys.size());
    }

    void clear_key_whitelist() {
        m_key_whitelist.clear();
        LOG_DEBUG("Cleared WriterRegistry key whitelist");
    }

private:
    std::unordered_map<std::type_index, WriterFn> m_writers;
    std::unordered_set<std::string> m_key_whitelist;
};

// inline WriterRegistry& global_registry() {
//   static WriterRegistry reg;
//   return reg;
// }

// template <class T>
// struct AutoRegister {
//   AutoRegister() {
//     global_registry().register_struct<T>();
//   }
// };