function gas = setsupersonicexitconditions(gas, expansionRatio)
import Gas.*
s_throat = State(gas);
s_exit = State();
for i = 1 : 1
    temperature_min = 0;
    temperature_max = s_throat.temperature;
    for j = 1 : 200
        setstate(gas, s_throat);
        s_exit.temperature = (temperature_min + temperature_max) / 2; % Bisection midpoint
        % Cantera hasn't implemented a function to change temperature at
        % constant entropy. I made my own function to faciliate that.
        settemperatureisentropic(gas, s_exit.temperature);
        s_exit = State(gas);

        massFlowFlux_calc = expansionRatio * s_exit.massFlowFlux;
        error_M = (massFlowFlux_calc - s_throat.massFlowFlux) / massFlowFlux_calc;
        if abs(error_M) < 1e-6  % Checking for convergence
            break;
        end
        if error_M < 0 % Range redefinition
            temperature_min = s_exit.temperature;
        else
            temperature_max = s_exit.temperature;
        end
    end
end
end