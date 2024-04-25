function solverOutputs = optimalmixture(solverInputs)

import Gas.*
import State.*

propellantArray = solverInputs.propellantArray;
gas = solverInputs.gas;
atmo = solverInputs.atmo;
expansionRatio = solverInputs.expansionRatio;
contractionRatio = solverInputs.contractionRatio;
heatMass = solverInputs.heatMass;
temperature_chamber_max = solverInputs.temperature_chamber_max;
mass_mole = solverInputs.mass_mole;

tolerance = 1e-9;
iterationLimit = 100;


% Point 1
ratio(1) = 1;
gas = combine(gas, propellantArray, [1, ratio(1)], mass_mole);
s_start = State(gas);
performance(1) = calculateperformance_optimal(gas, s_start, expansionRatio, contractionRatio, heatMass, atmo);

% Point 2
ratio(2) = 4;
gas = combine(gas, propellantArray, [1, ratio(2)], mass_mole);
s_start = State(gas);
performance(2) = calculateperformance_optimal(gas, s_start, expansionRatio, contractionRatio, heatMass, atmo);

% Point 3
ratio(3) = 7;
gas = combine(gas, propellantArray, [1, ratio(3)], mass_mole);
s_start = State(gas);
performance(3) = calculateperformance_optimal(gas, s_start, expansionRatio, contractionRatio, heatMass, atmo);

numericalMethod = LocalMinMax("optimalmixture_optimal", tolerance, 100, ratio, performance, 'max');
numericalMethod.setX_min_max(0, 10);
while 1
    ratio = numericalMethod.findnewX;
    gas = combine(gas, propellantArray, [1, ratio], mass_mole);
    s_start = State(gas);
    [performance, solverOutputs] = calculateperformance_optimal(gas, s_start, expansionRatio, contractionRatio, heatMass, atmo);

    numericalMethod.updateXY(ratio, performance);
    if numericalMethod.checkconvergence
        break;
    end
end
if solverOutputs.s_chamber.temperature > temperature_chamber_max
    % Point 1
    ratio(1) = 1;
    gas = combine(gas, propellantArray, [1, ratio(1)], mass_mole);
    s_start = State(gas);
    [~, ~, error_T(1)] = calculatetemperature(gas, s_start, expansionRatio, contractionRatio, heatMass, atmo, temperature_chamber_max);

    % Point 2
    ratio(2) = 4;
    gas = combine(gas, propellantArray, [1, ratio(2)], mass_mole);
    s_start = State(gas);
    [~, ~, error_T(2)] = calculatetemperature(gas, s_start, expansionRatio, contractionRatio, heatMass, atmo, temperature_chamber_max);

    % Point 3
    ratio(3) = 7;
    gas = combine(gas, propellantArray, [1, ratio(3)], mass_mole);
    s_start = State(gas);
    [~, ~, error_T(3)] = calculatetemperature(gas, s_start, expansionRatio, contractionRatio, heatMass, atmo, temperature_chamber_max);

    numericalMethod = MullersMethod("optimalmixture_maxtemp", tolerance, 100, ratio, error_T, 'min');
    numericalMethod.setX_min_max(0, 100);
    while 1
        ratio = numericalMethod.findnewX;
        gas = combine(gas, propellantArray, [1, ratio], mass_mole);
        s_start = State(gas);
        [performance, solverOutputs, error_T] = calculatetemperature(gas, s_start, expansionRatio, contractionRatio, heatMass, atmo, temperature_chamber_max);

        numericalMethod.updateXY(ratio, error_T);
        if numericalMethod.checkconvergence
            break;
        end
    end
end
% log = zeros(100,2);
% i = 0;
% for ratio = 2 : 0.05 : 6
%     i = i + 1;
%     gas = combine(gas, propellantArray, [1, ratio(1)], mass_mole);
%     s_start = State(gas);
%     [performance, solverOutputs] = calculateperformance_optimal(gas, s_start, expansionRatio, contractionRatio, heatMass, atmo);
%     log(i,:) = [ratio, performance];
% end 
solverOutputs.ratio = [1, ratio];
end

function [performance, solverOutputs, temperature] = calculateperformance_optimal(gas, s_injection, expansionRatio, contractionRatio, heatMass, atmo)
import Gas.*
setchamberconditions(gas, heatMass);
setarearatioisentropic(gas, contractionRatio, "sub");
s_chamber = State(gas);
s_throat = gas.sonic;
s_supersonicExit = setsupersonicexitconditions(gas, expansionRatio);

effectiveExhaustVelocity = effectiveexhaustvelocity(gas, atmo);
temperature = s_chamber.temperature;
performance = effectiveExhaustVelocity;
solverOutputs.s_chamber = s_chamber;
solverOutputs.s_throat = s_throat;
solverOutputs.s_supersonicExit = s_supersonicExit;
solverOutputs.effectiveExhaustVelocity = effectiveExhaustVelocity;
end

function [performance, solverOutputs, temperature] = calculatetemperature(gas, s_injection, expansionRatio, contractionRatio, heatMass, atmo, max_temp)
import Gas.*
setchamberconditions(gas, heatMass);
setarearatioisentropic(gas, contractionRatio, "sub");
s_chamber = State(gas);
s_throat = gas.sonic;
s_supersonicExit = setsupersonicexitconditions(gas, expansionRatio);
effectiveExhaustVelocity = effectiveexhaustvelocity(gas, atmo);
temperature = (max_temp - s_chamber.temperature) / max_temp;
performance = effectiveExhaustVelocity;
solverOutputs.s_chamber = s_chamber;
solverOutputs.s_throat = s_throat;
solverOutputs.s_supersonicExit = s_supersonicExit;
solverOutputs.effectiveExhaustVelocity = effectiveExhaustVelocity;
end