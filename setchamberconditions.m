% Author: Bartosz Kruszona, 2023

% Combustion
% Cantera solves for chemical equilibrium
% For non-chemical engines, the heat is added to the gas at constant
% pressure
function setchamberconditions(gas, s_injection, heatPower)
setstate(gas, s_injection);
import State.*
import Gas.*
equilibrate(gas.solution, 'HP');
setstate(gas, 'velocity', s_injection.massFlowFlux / gas.density);
if nargin == 2 || heatPower == 0
    return;
end

% Declaration of states.
s_preHeat = State(gas);
s_postHeat = State(gas);

% Error declarations.
Tolerance = 1e-12;
error_E = Tolerance * 10;
error_E_old = error_E;
error_E_abs = 0;

% Estimation of final temperature
heat_delta = heatPower / gas.massFlowFlux;

n = 0;
while abs(error_E) > Tolerance
    if abs(error_E_old / error_E) < 1
        n = n + 1;
        if n > 2
            disp("convergence failed in setchamberconditions with error = " + error_E);
            break;
        end
    else
        n = 0;
    end
    error_E_old = error_E;
    
    error_E_abs = error_E_abs + (s_preHeat.totEnergy + heat_delta - s_postHeat.totEnergy);
    s_postHeat.temperature = s_preHeat.temperature + (heat_delta + error_E_abs) / s_preHeat.cp;
    setstate(gas, 'T', s_postHeat.temperature, 'P', s_preHeat.pressure);
    setstate(gas, 'velocity', s_injection.massFlowFlux / gas.density);
    s_postHeat = State(gas);

    error_E = (s_preHeat.totEnergy + heat_delta - s_postHeat.totEnergy) / (s_preHeat.totEnergy + heat_delta); % Checking for convergence
end
end