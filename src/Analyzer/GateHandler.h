#ifndef GATE_HANDLER_H
#define GATE_HANDLER_H

#include <string>
#include <unordered_map>
#include <memory>

class GateHandler
{
public:
    using Ref = std::shared_ptr<GateHandler>;

    struct Gate
    {
        std::string name = "";
        double minValue = 0.0;
        double maxValue = 0.0;
    };

    GateHandler() {}
    ~GateHandler() {}

    void AddGate(const Gate& gate)
    {
        m_gateMap[gate.name] = gate;
    }

    bool IsInside(const std::string& name, double value)
    {
        auto iter = m_gateMap.find(name);
        if(iter != m_gateMap.end())
        {
            return value > iter->second.minValue && value < iter->second.maxValue;
        }
        return false;
    }

private:
    std::unordered_map<std::string, Gate> m_gateMap;
};

#endif