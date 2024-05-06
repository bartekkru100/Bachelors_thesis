function effectiveexhaustvelocity = effectiveexhaustvelocity(gas, s_atmo)

% Calculates the specific impulse
effectiveexhaustvelocity = (gas.velocity + (gas.pressure - s_atmo.pressure) / gas.massFlowFlux);
end