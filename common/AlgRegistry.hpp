#pragma once

#include <yaml-cpp/yaml.h>

#include <functional>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <unordered_map>

#include "common/IAlg.hpp"
#include "common/RunContext.hpp"

class AlgRegistry {
public:
    using Creator = std::function<std::unique_ptr<IAlg>(RunContext&, const YAML::Node&)>;

    static AlgRegistry& instance() {
        static AlgRegistry r;
        return r;
    }

    void add(std::string type, Creator c) {
        std::lock_guard<std::mutex> lock(m_);
        auto [it, ok] = creators_.emplace(std::move(type), std::move(c));
        if (!ok) {
            LOG_ERROR("AlgRegistry: duplicate registration for type: {}", it->first);
            throw std::runtime_error("AlgRegistry: duplicate registration for type: " + it->first);
        }
    }

    std::unique_ptr<IAlg> create(const std::string& type, RunContext& ctx, const YAML::Node& cfg) const {
        std::lock_guard<std::mutex> lock(m_);
        auto it = creators_.find(type);
        if (it == creators_.end()) {
            LOG_ERROR("AlgRegistry: unknown alg type: {}", type);
            throw std::runtime_error("AlgRegistry: unknown alg type: " + type);
        }
        return (it->second)(ctx, cfg);
    }

private:
    AlgRegistry() = default;

    mutable std::mutex m_;
    std::unordered_map<std::string, Creator> creators_;
};

namespace algreg_detail {

template <class T>
struct Registrar {
    explicit Registrar(const char* type_name, std::function<void(T&)> post_init = nullptr) {
        AlgRegistry::instance().add(
            type_name,
            [type_name, post_init](RunContext& ctx, const YAML::Node& cfg) -> std::unique_ptr<IAlg> {
                auto alg = std::make_unique<T>(ctx, std::string(type_name));
                alg->parse_cfg(cfg);
                if (post_init) post_init(*alg);
                return alg;
            });
    }
};

} // namespace algreg_detail

// helper for unique static variable names (Typeに :: があっても壊れない)
#define AHCAL_CONCAT_IMPL(a, b) a##b
#define AHCAL_CONCAT(a, b) AHCAL_CONCAT_IMPL(a, b)

#define AHCAL_REGISTER_ALG(Type, NameLiteral) namespace { const ::algreg_detail::Registrar<Type> AHCAL_CONCAT(kAlgReg_, __COUNTER__)(NameLiteral); }

#define AHCAL_REGISTER_ALG_POSTINIT(Type, NameLiteral, PostInitLambda) namespace { const ::algreg_detail::Registrar<Type> AHCAL_CONCAT(kAlgReg_, __COUNTER__)(NameLiteral, PostInitLambda); }