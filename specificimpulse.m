function specificImpulse = specificimpulse(gas, s_atmo, isInSeconds)

% Calculates the specific impulse
if isInSeconds == false
specificImpulse = (gas.velocity + (gas.pressure - s_atmo.pressure) / gas.massFlowFlux);
else
specificImpulse = (gas.velocity + (gas.pressure - s_atmo.pressure) / gas.massFlowFlux) / 9.81;
end
end