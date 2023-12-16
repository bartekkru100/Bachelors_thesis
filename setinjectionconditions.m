function gas = setinjectionconditions(gas, s_injection, contractionRatio)
import Gas.*
massFlowFlux_temp = gas.massFlowFlux;
setstate(gas, s_injection);
setstate(gas, 'velocity', massFlowFlux_temp / (gas.density * contractionRatio) );
end