function specificImpulse = specificimpulse(gas, s_atmo)
specificImpulse = (gas.velocity + (gas.pressure - s_atmo.pressure) / gas.massFlowFlux);% / 9.81
end