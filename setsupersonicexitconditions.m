function [expansionRatio, error_M] = setsupersonicexitconditions(gas, s_throat, expansionRatio)
import Gas.*

% This function calculates the exit state assuming a completely supersonic
% flow in the diverging part of the nozzle. Iterating over pressure to
% satisfy mass conservation.

% Using Muller's method

tolerance = 1e-12;
iterationLimit = 100;

% Point 1
pressure(1) = s_throat.pressure / expansionRatio ^ 1.4 / 2;
error_M(1) = error_exit(gas, pressure(1), s_throat, expansionRatio);

% Point 2
pressure(2) = s_throat.pressure / expansionRatio ^ 1.2 / 1.5;
error_M(2) = error_exit(gas, pressure(2), s_throat, expansionRatio);

% Point 3
pressure(3) = s_throat.pressure / expansionRatio;
error_M(3) = error_exit(gas, pressure(3), s_throat, expansionRatio);

numericalMethod = MullersMethod("setsupersonicexitconditions", tolerance, iterationLimit, pressure, error_M, 'min');
numericalMethod.setX_min(0);

while 1
    pressure = numericalMethod.findnewX;

    error_M = error_exit(gas, pressure, s_throat, expansionRatio);
    
    numericalMethod.updateXY(pressure, error_M);
    if numericalMethod.checkconvergence
        break;
    end
end
end

function error_M = error_exit(gas, pressure, s_throat, expansionRatio)
import Gas.*

setstate(gas, s_throat);
setpressureisentropic(gas, pressure, s_throat);
s_supersonicExit = State(gas);

massFlowFlux_calc = expansionRatio * s_supersonicExit.massFlowFlux;
error_M = (massFlowFlux_calc - s_throat.massFlowFlux) / massFlowFlux_calc;
end