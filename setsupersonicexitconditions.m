function gas = setsupersonicexitconditions(gas, expansionRatio)
import Gas.*

s_throat = State(gas);
Tolerance = 1e-12;

% Using Muller's method

% Point 1
temperature_1 = s_throat.temperature / sqrt(expansionRatio);
error_M_1 = exiterror(gas, temperature_1, s_throat, expansionRatio);

% Point 2
temperature_2 = s_throat.temperature / 2;
error_M_2 = exiterror(gas, temperature_2, s_throat, expansionRatio);

% Point 3
temperature_3 = s_throat.temperature;
error_M_3 = exiterror(gas, temperature_3, s_throat, expansionRatio);

n = 0;
while 1
    n = n + 1;
    if n > 100
        disp("convergence failed")
        break;
    end
    q = (temperature_1 - temperature_2) / (temperature_2 - temperature_3);
    a = q * error_M_1 - q * (1 + q) * error_M_2 + q ^ 2 * error_M_3;
    b = (2 * q + 1) * error_M_1 - (1 + q) ^ 2 * error_M_2 + q ^ 2 * error_M_3;
    c = (1 + q) * error_M_1;
    sqrtDelta = sqrt(b ^ 2 - 4 * a * c);
    denom(1) = (b + sqrtDelta);
    denom(2) = (b - sqrtDelta);
    temperature_0 = min(temperature_1 - (temperature_1 - temperature_2) * (2 * c) ./ denom);

    if temperature_0 <= 0
        temperature_0 = min([temperature_1, temperature_2, temperature_3]) / 3;
        if temperature_1 <= 0
            disp("error")
            break;
        end
    end

    error_M_0 = exiterror(gas, temperature_0, s_throat, expansionRatio);

    temperature_3 = temperature_2;
    error_M_3 = error_M_2;
    temperature_2 = temperature_1;
    error_M_2 = error_M_1;
    temperature_1 = temperature_0;
    error_M_1 = error_M_0;

    if abs(error_M_1) < Tolerance
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