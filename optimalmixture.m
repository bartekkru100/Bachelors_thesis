function solverOutputs = optimalmixture(solverInputs)

import Gas.*
import State.*

propellantArray = solverInputs.propellantArray;
atmo = solverInputs.atmo;
expansionRatio = solverInputs.expansionRatio;
contractionRatio = solverInputs.contractionRatio;
heatPowerArea = solverInputs.heatPowerArea;
temperature_chamber_max = solverInputs.temperature_chamber_max;

tolerance = 1e-9;
iterationLimit = 100;

mode = 'mass';

% Point 1
ratio(1) = 1;
gas = combine('new', propellantArray, [1, ratio(1)], mode);
s_start = State(gas);
value(1) = calculateperformance_optimal(gas, s_start, expansionRatio, heatPowerArea, atmo);

% Point 2
ratio(2) = 4;
gas = combine(gas, propellantArray, [1, ratio(2)], mode);
s_start = State(gas);
value(2) = calculateperformance_optimal(gas, s_start, expansionRatio, heatPowerArea, atmo);

% Point 3
ratio(3) = 7;
gas = combine(gas, propellantArray, [1, ratio(3)], mode);
s_start = State(gas);
value(3) = calculateperformance_optimal(gas, s_start, expansionRatio, heatPowerArea, atmo);

numericalMethod = LocalMinMax("optimalmixture_optimal", tolerance, 100, ratio, value, 'max');
numericalMethod.setX_min_max(0, 10);
while 1
    ratio = numericalMethod.findnewX;
    gas = combine(gas, propellantArray, [1, ratio], mode);
    s_start = State(gas);
    [value, s_chamber] = calculateperformance_optimal(gas, s_start, expansionRatio, heatPowerArea, atmo);

    numericalMethod.updateXY(ratio, value);
    if numericalMethod.checkconvergence
        break;
    end
end
if s_chamber.temperature > temperature_chamber_max
    % Point 1
    ratio(1) = 1;
    gas = combine(gas, propellantArray, [1, ratio(1)], mode);
    s_start = State(gas);
    error_T(1) = calculatetemperature(gas, s_start, expansionRatio, contractionRatio, heatPowerArea, atmo, temperature_chamber_max);

    % Point 2
    ratio(2) = 4;
    gas = combine(gas, propellantArray, [1, ratio(2)], mode);
    s_start = State(gas);
    error_T(2) = calculatetemperature(gas, s_start, expansionRatio, contractionRatio, heatPowerArea, atmo, temperature_chamber_max);

    % Point 3
    ratio(3) = 7;
    gas = combine(gas, propellantArray, [1, ratio(3)], mode);
    s_start = State(gas);
    error_T(3) = calculatetemperature(gas, s_start, expansionRatio, contractionRatio, heatPowerArea, atmo, temperature_chamber_max);

    numericalMethod = MullersMethod("optimalmixture_maxtemp", tolerance, 100, ratio, error_T, 'min');
    numericalMethod.setX_min_max(0, 10);
    while 1
        ratio = numericalMethod.findnewX;
        gas = combine(gas, propellantArray, [1, ratio], mode);
        s_start = State(gas);
        error_T = calculatetemperature(gas, s_start, expansionRatio, contractionRatio, heatPowerArea, atmo, temperature_chamber_max);

        numericalMethod.updateXY(ratio, error_T);
        if numericalMethod.checkconvergence
            break;
        end
    end
end
ratio

i = 1;
for ratio = 1:0.05:6
    i = i + 1;
    gas = combine(gas, propellantArray, [1, ratio], mode);
    s_start = State(gas);
    log(i, :) = [ratio, calculateperformance_optimal(gas, s_start, expansionRatio, heatPowerArea, atmo)];
end


end

function [specificImpulse_ms, s_chamber] = calculateperformance_optimal(gas, s_start, expansionRatio, heatPowerArea, atmo)
tolerance = 1e-9;
s_injection = s_start;
setchamberconditions(gas, s_injection)
s_chamber = State(gas);
[s_throat, s_stagnation] = gas.sonic;
s_supersonicExit = setsupersonicexitconditions(gas, s_throat, expansionRatio);
specificImpulse_ms = specificimpulse(gas, atmo, false);
velocity_star = s_chamber.pressure / s_throat.massFlowFlux;
%output = s_supersonicExit.velocity;
output = specificImpulse_ms;
end

function temperature = calculatetemperature(gas, s_start, expansionRatio, contractionRatio, heatPowerArea, atmo, max_temp)
tolerance = 1e-9;
s_injection = s_start;
setchamberconditions(gas, s_injection);
s_chamber = State(gas);
[s_throat, s_stagnation] = gas.sonic;
s_supersonicExit = setsupersonicexitconditions(gas, s_throat, expansionRatio);
temperature = (max_temp - s_chamber.temperature) / max_temp;
end