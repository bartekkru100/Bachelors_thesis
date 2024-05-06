function thrustCoefficient = thrustcoefficient(gas, s_atmo)
    thrustCoefficient = effectiveexhaustvelocity(gas, s_atmo) / characteristicvelocity(gas);
end