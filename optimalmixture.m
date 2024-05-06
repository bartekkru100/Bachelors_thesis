function solverOutputs = optimalmixture(solverInputs)

import Gas.*
import State.*

propellantArray = solverInputs.propellantArray;
gas = solverInputs.gas;
s_injection.pressure = solverInputs.pressure_injection;
s_injection.temperature = solverInputs.temperature_injection;
ratios = solverInputs.ratios;
atmo = solverInputs.atmo;
expansionRatio = solverInputs.expansionRatio;
contractionRatio = solverInputs.contractionRatio;
heatMass = solverInputs.heatMass;
thrust = solverInputs.thrust;
temperature_chamber_max = solverInputs.temperature_chamber_max;

mass_mole = solverInputs.mass_mole;
searchedIndicator = solverInputs.searchedIndicator;

tolerance = 1e-9;
iterationLimit = 100;

overheatingFound = false;
flowState = "Optimal mixture found";

nProps = length(ratios);
otherRatios = ratios(2:end);
otherRatios = otherRatios / sum(otherRatios);

if mass_mole == "mass"
    molecularWeightAdjustment = sqrt(propellantArray(1).meanMolecularWeight / propellantArray(2).meanMolecularWeight);
else
    molecularWeightAdjustment = 1;
end

% Point 1
ratio(1) = 0.25 * molecularWeightAdjustment;
setstate(gas, 'P', s_injection.pressure, 'T', s_injection.temperature);
gas = combine(gas, propellantArray, [ratio(1), otherRatios], mass_mole);
s_start = State(gas);
performance(1) = calculateperformance_optimal(gas, s_start, expansionRatio, contractionRatio, heatMass, atmo, searchedIndicator);

% Point 2
ratio(2) = 1 * molecularWeightAdjustment;
setstate(gas, 'P', s_injection.pressure, 'T', s_injection.temperature);
gas = combine(gas, propellantArray, [ratio(2), otherRatios], mass_mole);
s_start = State(gas);
performance(2) = calculateperformance_optimal(gas, s_start, expansionRatio, contractionRatio, heatMass, atmo, searchedIndicator);

% Point 3
ratio(3) = 4 * molecularWeightAdjustment;
setstate(gas, 'P', s_injection.pressure, 'T', s_injection.temperature);
gas = combine(gas, propellantArray, [ratio(3), otherRatios], mass_mole);
s_start = State(gas);
performance(3) = calculateperformance_optimal(gas, s_start, expansionRatio, contractionRatio, heatMass, atmo, searchedIndicator);

numericalMethod = LocalMinMax("optimalmixture_optimal", tolerance, 100, ratio, performance, 'max');
numericalMethod.setX_min(0);
while 1
    ratio = numericalMethod.findnewX;
    setstate(gas, 'P', s_injection.pressure, 'T', s_injection.temperature);
    gas = combine(gas, propellantArray, [ratio, otherRatios], mass_mole);
    s_start = State(gas);
    [performance, solverOutputs] = calculateperformance_optimal(gas, s_start, expansionRatio, contractionRatio, heatMass, atmo , searchedIndicator);

    numericalMethod.updateXY(ratio, performance);
    if numericalMethod.checkconvergence
        break;
    end
end
if solverOutputs.s_chamber.temperature > temperature_chamber_max
    overheatingFound = true;
    flowState = "Optimal mixture changed due to overheating";

    % Point 1
    ratio(1) = 0.5 * ratio;
    setstate(gas, 'P', s_injection.pressure, 'T', s_injection.temperature);
    gas = combine(gas, propellantArray, [ratio, otherRatios], mass_mole);
    s_start = State(gas);
    [~, ~, error_T(1)] = calculatetemperature(gas, s_start, expansionRatio, contractionRatio, heatMass, atmo, temperature_chamber_max, searchedIndicator);

    % Point 2
    ratio(2) = ratio;
    setstate(gas, 'P', s_injection.pressure, 'T', s_injection.temperature);
    gas = combine(gas, propellantArray, [ratio(2), otherRatios], mass_mole);
    s_start = State(gas);
    [~, ~, error_T(2)] = calculatetemperature(gas, s_start, expansionRatio, contractionRatio, heatMass, atmo, temperature_chamber_max, searchedIndicator);

    % Point 3
    ratio(3) = 2 * ratio;
    setstate(gas, 'P', s_injection.pressure, 'T', s_injection.temperature);
    gas = combine(gas, propellantArray, [ratio(3), otherRatios], mass_mole);
    s_start = State(gas);
    [~, ~, error_T(3)] = calculatetemperature(gas, s_start, expansionRatio, contractionRatio, heatMass, atmo, temperature_chamber_max, searchedIndicator);

    numericalMethod = MullersMethod("optimalmixture_maxtemp", tolerance, 100, ratio, error_T, 'min');
    numericalMethod.setX_min_max(0, 100);
    while 1
        ratio = numericalMethod.findnewX;
        setstate(gas, 'P', s_injection.pressure, 'T', s_injection.temperature);
        gas = combine(gas, propellantArray, [ratio, otherRatios], mass_mole);
        s_start = State(gas);
        [performance, solverOutputs, error_T] = calculatetemperature(gas, s_start, expansionRatio, contractionRatio, heatMass, atmo, temperature_chamber_max, searchedIndicator);

        numericalMethod.updateXY(ratio, error_T);
        if numericalMethod.checkconvergence
            break;
        end
    end
end
if numericalMethod.hasFailed
    global solverError;
    solverError = "ERROR: Convergence failed, try changing the pressure ratio.";
    error("Convergence failed, try changing the pressure ratio.")
end

s_injection = solverOutputs.s_injection;
s_chamber = solverOutputs.s_chamber;
s_throat = solverOutputs.s_throat;
s_exit = solverOutputs.s_exit;

effectiveExhaustVelocity = effectiveexhaustvelocity(gas, atmo);
characteristicVelocity = characteristicvelocity(gas, atmo);
thrustCoefficient = thrustcoefficient(gas, atmo);
specificImpulse = specificimpulse(gas, atmo);

solverOutputs.effectiveExhaustVelocity = effectiveExhaustVelocity;
solverOutputs.specificImpulse = specificImpulse;
solverOutputs.thrustCoefficient = thrustCoefficient;
solverOutputs.characteristicVelocity = characteristicVelocity;

massFlow = thrust / effectiveExhaustVelocity;
A_chamber = massFlow / s_chamber.massFlowFlux;
A_throat = massFlow / s_throat.massFlowFlux;
A_exit = massFlow / s_exit.massFlowFlux;
heatPower = heatMass * massFlow;
D_chamber = sqrt(4 * A_chamber / pi);
D_throat = sqrt(4 * A_throat / pi);
D_exit = sqrt(4 * A_exit / pi);

solverOutputs.A_chamber = A_chamber;
solverOutputs.A_throat = A_throat;
solverOutputs.A_exit = A_exit;
solverOutputs.D_chamber = D_chamber;
solverOutputs.D_throat = D_throat;
solverOutputs.D_exit = D_exit;
solverOutputs.massFlow = massFlow;
solverOutputs.heatPower = heatPower;

solverOutputs.shockPosition = [];
solverOutputs.s_shock_1 = [];
solverOutputs.s_shock_2 = [];

solverOutputs.s_injection = s_injection;
solverOutputs.s_chamber = s_chamber;
solverOutputs.s_throat = s_throat;
solverOutputs.s_exit = s_exit;
solverOutputs.optimalRatio = ratio;
solverOutputs.pressure_ideal = [];
solverOutputs.pressure_shockAtExit = [];
solverOutputs.pressure_shockAtThroat = [];
solverOutputs.pressure_separation = [];
solverOutputs.hasCondensation = gas.hasCondensation;
solverOutputs.expansionRatio_shockAtThroat = [];
solverOutputs.expansionRatio_optimal = [];
solverOutputs.expansionRatio_shockAtExit_1 = [];
solverOutputs.expansionRatio_shockAtExit_2 = [];
solverOutputs.expansionRatio_condensation = [];
solverOutputs.expansionRatio_separation = [];
solverOutputs.flowState = flowState;
end

function [performance, solverOutputs, temperature] = calculateperformance_optimal(gas, s_injection, expansionRatio, contractionRatio, heatMass, atmo, searchedIndicator)
import Gas.*
setchamberconditions(gas, heatMass);
s_chamber = State(gas);
s_throat = gas.sonic;
s_exit = setsupersonicexitconditions(gas, expansionRatio);

if searchedIndicator == "specificImpulse"
    performance = effectiveexhaustvelocity(gas, atmo);
elseif searchedIndicator == "characteristicVelocity"
    performance = characteristicvelocity(gas);
else
    performance = thrustcoefficient(gas, atmo);
end
temperature = s_chamber.temperature;
solverOutputs.s_injection = s_injection;
solverOutputs.s_chamber = s_chamber;
solverOutputs.s_throat = s_throat;
solverOutputs.s_exit = s_exit;
end

function [performance, solverOutputs, temperature] = calculatetemperature(gas, s_injection, expansionRatio, contractionRatio, heatMass, atmo, max_temp, searchedIndicator)
import Gas.*
setchamberconditions(gas, heatMass);
s_chamber = State(gas);
s_throat = gas.sonic;
s_exit = setsupersonicexitconditions(gas, expansionRatio);

if searchedIndicator == "specificImpulse"
    performance = effectiveexhaustvelocity(gas, atmo);
elseif searchedIndicator == "characteristicVelocity"
    performance = characteristicvelocity(gas);
else
    performance = thrustcoefficient(gas, atmo);
end
temperature = (max_temp - s_chamber.temperature) / max_temp;
solverOutputs.s_injection = s_injection;
solverOutputs.s_chamber = s_chamber;
solverOutputs.s_throat = s_throat;
solverOutputs.s_exit = s_exit;
end