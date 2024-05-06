function characteristicVelocity = characteristicvelocity(gas, s_atmo)
    characteristicVelocity = gas.stagnation.pressure / gas.sonic.massFlowFlux;
end