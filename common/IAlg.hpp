#pragma once
#include "common/EventStore.hpp"
#include <yaml-cpp/yaml.h>
struct RunContext;

class IAlg {
public:
    IAlg(RunContext& ctx, std::string name) : m_ctx(ctx), m_name(std::move(name)) {}
    virtual ~IAlg() = default;

    virtual void initialize() {}   // optional
    virtual void execute(EventStore& evt) = 0;
    virtual void finalize() {}     // optional
    virtual void  parse_cfg(const YAML::Node& n) = 0;
protected:
    RunContext& ctx() { return m_ctx; }
    const RunContext& ctx() const { return m_ctx; }
    const std::string& name() const { return m_name; }

private:
    RunContext& m_ctx;
    std::string m_name;
};
