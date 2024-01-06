function gas = setsupersonicexitconditions(gas, expansionRatio)
import Gas.*

s_throat = State(gas);
tolerance = 1e-12;

% Using Muller's method

% Point 1
temperature(1) = s_throat.temperature / sqrt(expansionRatio);
error_M(1) = exiterror(gas, temperature(1), s_throat, expansionRatio);

% Point 2
temperature(2) = s_throat.temperature / 2;
error_M(2) = exiterror(gas, temperature(2), s_throat, expansionRatio);

% Point 3
temperature(3) = s_throat.temperature;
error_M(3) = exiterror(gas, temperature(3), s_throat, expansionRatio);

numericalMethod = MullersMethod("setsupersonicexitconditions", tolerance, 100, temperature, error_M, 'min');

while 1
    temperature_0 = numericalMethod.findnewX;

    if temperature_0 <= 0
        temperature_0 = min([temperature(1), temperature(2), temperature(3)]) / 3;
        if temperature(1) <= 0
            disp("error")
            break;
        end
    end

    error_M_0 = exiterror(gas, temperature_0, s_throat, expansionRatio);
    
    numericalMethod.updateXY(temperature_0, error_M_0);
    if numericalMethod.checkconvergence
        break;
    end
end
end

function error_M = exiterror(gas, temperature_2, s_throat, expansionRatio)
import Gas.*
setstate(gas, s_throat);
settemperatureisentropic(gas, temperature_2);
s_exit = State(gas);

massFlowFlux_calc = expansionRatio * s_exit.massFlowFlux;
error_M = (massFlowFlux_calc - s_throat.massFlowFlux) / massFlowFlux_calc;
end