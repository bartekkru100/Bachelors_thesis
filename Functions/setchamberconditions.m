function setchamberconditions(gas, s_injection, heatPower)
import State.*
import Gas.*

% This sets the chemical equilibrium and apllies heat from a heating
% element

if ~(nargin == 2 || heatPower == 0)

    s_preHeat = State(gas);

    % Error declarations.
    tolerance = 1e-9;
    maxIterations = 10;
    error_E_abs = 0;

    heat_delta = heatPower / gas.massFlowFlux;

    numericalMethod = ErrorCorrectedMethod("setchamberconditions", tolerance, maxIterations);
    numericalMethod.setX_min(0);

    while 1
        error_E_abs = error_E_abs + (s_preHeat.totEnergy + heat_delta - gas.totEnergy);
        temperature = s_preHeat.temperature + (heat_delta + error_E_abs) / s_preHeat.cp;

        temperature = numericalMethod.findnewX(temperature);

        error_E = error_heatingelement(gas, temperature, s_preHeat, s_injection, heat_delta);

        numericalMethod.updateXY(temperature, error_E);
        if numericalMethod.checkconvergence
            break;
        end
    end
end
equilibrate(gas.solution, 'HP');
%setisentropic
setstate(gas, 'velocity', s_injection.massFlowFlux / gas.density);
end

%--------------------------------------------------------------------------

% The function checks for energy conservation

function error_E = error_heatingelement(gas, temperature, s_preHeat, s_injection, heat_delta)
import Gas.*

setstate(gas, 'T', temperature, 'P', s_preHeat.pressure);
setstate(gas, 'velocity', s_injection.massFlowFlux / gas.density);

error_E = (s_preHeat.totEnergy + heat_delta - gas.totEnergy) / (s_preHeat.totEnergy + heat_delta);
end