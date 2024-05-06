function specificImpulse = specificimpulse(gas, s_atmo)

% Calculates the specific impulse
specificImpulse = (gas.velocity + (gas.pressure - s_atmo.pressure) / gas.massFlowFlux) / 9.80665;
end