import Gas.*
import Unit.*
import PhaseDiagram.*
import State.*
path_Classes = [pwd, '\Classes'];
path_Functions = [pwd, '\Functions'];

addpath(path_Classes, path_Functions)
% Test file to test parts of code


clc, clear, cleanup

atmo = Gas('air');
s_atmo = State(atmo);
s_atmo.pressure = 1;
setstate(atmo, s_atmo);
s_atmo = State(atmo);
expansionRatio = 77.5;
contractionRatio = 15;
tolerance = 1e-9;
iterationLimit = 100;
heatPowerArea = 0;
thrust = 1.86e6;
efficiency = 1;
temperature_max = 3500;

tic

%setstate(gas, 'Y', 'H2:1,O2:6.03', 'T', 150, 'P', 20.64e6);
gas1 = Gas('R1highT');
setstate(gas1, 'Y', 'O2:1', 'T', 150, 'P', 20.64e6);
gas2 = Gas('R1highT');
setstate(gas2, 'Y', 'H2:1', 'T', 150, 'P', 20.64e6);

% Point 1
ratio(1) = 0.75;
gas = combine('new', [gas1, gas2], [ratio(1), 1 - ratio(1)]);
s_start = State(gas);
value(1) = calculateperformance_optimal(gas, s_start, expansionRatio, heatPowerArea, s_atmo);

% Point 2
ratio(2) = 0.5;
gas = combine(gas, [gas1, gas2], [ratio(2), 1 - ratio(2)]);
s_start = State(gas);
value(2) = calculateperformance_optimal(gas, s_start, expansionRatio, heatPowerArea, s_atmo);

% Point 3
ratio(3) = 0.25;
gas = combine(gas, [gas1, gas2], [ratio(3), 1 - ratio(3)]);
s_start = State(gas);
value(3) = calculateperformance_optimal(gas, s_start, expansionRatio, heatPowerArea, s_atmo);

numericalMethod = LocalMinMax("optimalmixture_optimal", tolerance, 100, ratio, value, 'max');
numericalMethod.setX_min_max(0, 1);
while 1
    ratio_0 = numericalMethod.findnewX;
    gas = combine(gas, [gas1, gas2], [ratio_0, 1 - ratio_0]);
    s_start = State(gas);
    value_0 = calculateperformance_optimal(gas, s_start, expansionRatio, heatPowerArea, s_atmo);

    numericalMethod.updateXY(ratio_0, value_0);
    if numericalMethod.checkconvergence
        break;
    end
end
toc
%{
tic
% Point 1
ratio(1) = 0.75;
gas = combine('new', [gas1, gas2], [ratio(1), 1 - ratio(1)]);
s_start = State(gas);
error_T(1) = calculatetemperature(gas, s_start, expansionRatio, contractionRatio, heatPowerArea, s_atmo, temperature_max);

% Point 2
ratio(2) = 0.5;
gas = combine(gas, [gas1, gas2], [ratio(2), 1 - ratio(2)]);
s_start = State(gas);
error_T(2) = calculatetemperature(gas, s_start, expansionRatio, contractionRatio, heatPowerArea, s_atmo, temperature_max);

% Point 3
ratio(3) = 0.25;
gas = combine(gas, [gas1, gas2], [ratio(3), 1 - ratio(3)]);
s_start = State(gas);
error_T(3) = calculatetemperature(gas, s_start, expansionRatio, contractionRatio, heatPowerArea, s_atmo, temperature_max);

numericalMethod = MullersMethod("optimalmixture_maxtemp", tolerance, 100, ratio, error_T, 'min');
numericalMethod.setX_min_max(0, 1);
while 1
    ratio_0 = numericalMethod.findnewX;
    gas = combine(gas, [gas1, gas2], [ratio_0, 1 - ratio_0]);
    s_start = State(gas);
    error_T_0 = calculatetemperature(gas, s_start, expansionRatio, contractionRatio, heatPowerArea, s_atmo, temperature_max);

    numericalMethod.updateXY(ratio_0, error_T_0);
    if numericalMethod.checkconvergence
        break;
    end
end
toc
%}
%{
i = 1;
for ratio = 0.2:0.05:0.5
    i = i + 1;
    gas = combine(gas, [gas1, gas2], [ratio, 1 - ratio]);
    s_start = State(gas);
    log(i, :) = [ratio, calculateperformance_optimal(gas, s_start, expansionRatio, heatPowerArea, s_atmo)]
end
%}



function specificImpulse_ms = calculateperformance_optimal(gas, s_start, expansionRatio, heatPowerArea, s_atmo)
tolerance = 1e-9;
s_injection = s_start;
setchamberconditions(gas, s_injection)
s_chamber = State(gas);
[s_throat, s_stagnation] = gas.sonic;
s_supersonicExit = setsupersonicexitconditions(gas, s_throat, expansionRatio);
specificImpulse_ms = specificimpulse(gas, s_atmo, false);
end

function temperature = calculatetemperature(gas, s_start, expansionRatio, contractionRatio, heatPowerArea, s_atmo, max_temp)
tolerance = 1e-9;
s_injection = s_start;
setchamberconditions(gas, s_injection)
s_chamber = State(gas);
[s_throat, s_stagnation] = gas.sonic;
s_supersonicExit = setsupersonicexitconditions(gas, s_throat, expansionRatio);
temperature = (max_temp - s_chamber.temperature) / max_temp;
end

function displaygasinfo(gas, name, areaRatio)
disp(name + ':' + ...
    " Velocity: " + gas.velocity + ...
    " Entropy: " + gas.entropy + ...
    " Total Energy: " + gas.totEnergy + ...
    " Mass Flow Rate: " + gas.massFlowFlux * areaRatio)
end
