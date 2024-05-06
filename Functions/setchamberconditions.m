function setchamberconditions(gas, contractionRatio, heatMass)
import State.*
import Gas.*

s_preHeat = State(gas);

% This sets the chemical equilibrium and apllies heat from a heating
% element

if ~(nargin == 1 || heatMass == 0)

    % Error declarations.
    tolerance = 1e-9;
    maxIterations = 10;
    error_E_abs = -heatMass;

    heat_delta = heatMass;

    gas.stagnation(false);
    setstate(gas, 'P', s_preHeat.pressure, 'H', s_preHeat.enthalpy + heat_delta);

    % numericalMethod = ErrorCorrectedMethod("setchamberconditions_stage_2", tolerance, maxIterations);
    % numericalMethod.disablewarnings;
    % numericalMethod.setX_min(0);
    % 
    % while 1
    %     error_E_abs = error_E_abs + (s_preHeat.totEnergy + heat_delta - gas.totEnergy);
    %     temperature = s_preHeat.temperature + (heat_delta + error_E_abs) / s_preHeat.cp;
    % 
    %     temperature = numericalMethod.findnewX(temperature);
    % 
    %     error_E = error_heatingelement(gas, temperature, s_preHeat, heat_delta);
    % 
    %     numericalMethod.updateXY(temperature, error_E);
    %     if numericalMethod.checkconvergence
    %         break;
    %     end
    % end
    % if numericalMethod.hasFailed
    % 
    %     tolerance = 1e-9;
    %     maxIterations = 10;
    % 
    %     temperature(1) = s_preHeat.temperature + (heat_delta) / s_preHeat.cp;
    %     error_E(1) = error_heatingelement(gas, temperature(1), s_preHeat, heat_delta);
    % 
    %     temperature(2) = s_preHeat.temperature + (heat_delta) / s_preHeat.cp / 2;
    %     error_E(2) = error_heatingelement(gas, temperature(2), s_preHeat, heat_delta);
    % 
    %     temperature(3) = s_preHeat.temperature + (heat_delta) / s_preHeat.cp / 3;
    %     error_E(3) = error_heatingelement(gas, temperature(3), s_preHeat, heat_delta);
    % 
    %     numericalMethod = MullersMethod("setchamberconditions_stage_2", tolerance, maxIterations, temperature, error_E, "max");
    %     numericalMethod.setX_min(0);
    % 
    %     while 1
    %         temperature = numericalMethod.findnewX;
    % 
    %         error_E = error_heatingelement(gas, temperature, s_preHeat, heat_delta);
    % 
    %         numericalMethod.updateXY(temperature, error_E);
    %         if numericalMethod.checkconvergence
    %             break;
    %         end
    %     end
    % end
end
s_postHeat = State(gas);
gas.equilibrium;
setarearatioisentropic(gas, contractionRatio, "sub");

end

%--------------------------------------------------------------------------

% The function checks for energy conservation

function error_E = error_heatingelement(gas, temperature, s_preHeat, heat_delta)
import Gas.*

setstate(gas, 'T', temperature, 'P', s_preHeat.pressure);
setstate(gas, 'velocity', s_preHeat.massFlowFlux / gas.density);

error_E = (s_preHeat.totEnergy + heat_delta - gas.totEnergy) / (s_preHeat.totEnergy + heat_delta);
end