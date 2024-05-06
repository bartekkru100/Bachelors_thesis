function solverOutputs = engineanalysis(solverInputs)

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

separationTolerance = solverInputs.separationTolerance;
mass_mole = solverInputs.mass_mole;

setstate(gas, 'P', s_injection.pressure, 'T', s_injection.temperature);
gas = combine(gas, propellantArray, ratios, mass_mole);
s_injection = State(gas);
setchamberconditions(gas, contractionRatio, heatMass);
s_chamber = State(gas);
s_throat = gas.sonic;
s_exit = setsupersonicexitconditions(gas, expansionRatio);
s_stagnation = gas.stagnation;

pressure_shockAtThroat = findshockatthroat_pressure(gas, expansionRatio);
pressure_shockAtExit = findshockatexit_pressure(gas, s_exit, expansionRatio);
pressure_separation  = findseparation_pressure(gas, s_exit, separationTolerance);
pressure_ideal = findidealexpansion_pressure(gas, s_exit);


shockFound = false;
separationFound = false;

if atmo.pressure > s_stagnation.pressure
    flowState = "Stagation pressure below ambient, flow impossible";

    s_injection = State(setstate(gas, 'P', atmo.pressure, 'T', atmo.temperature, 'velocity', 0));
    s_chamber = s_injection;
    s_throat = s_injection;
    s_exit = s_injection;
else
    if atmo.pressure > pressure_shockAtThroat
        [s_exit, s_throat, s_chamber] = setsubsonicexitconditions(gas, atmo, expansionRatio, contractionRatio);
        flowState.current = "Fully subsonic";
    else
        if atmo.pressure > pressure_shockAtExit
            [shockPosition, s_shock_1, s_shock_2, s_exit] = shockposition(gas, atmo, s_exit);
            flowState = "Normal shock inside the nozzle";
            shockFound = true;
        else
            if atmo.pressure >= pressure_ideal
                flowState = "Fully supersonic, overexpanded";
                if atmo.pressure > pressure_separation
                    separationFound = true;
                    flowState = flowState + ", separation likely";
                end
            else
                flowState = "Fully supersonic, underexpanded";
            end
        end
    end
end

if gas.hasCondensation
    flowState = flowState + ", condensation in the nozzle likely";
end

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
if s_chamber.temperature > solverInputs.temperature_chamber_max
    flowState = flowState + ", chamber overheating";
else
end
if shockFound
    solverOutputs.shockPosition = shockPosition;
    solverOutputs.s_shock_1 = s_shock_1;
    solverOutputs.s_shock_2 = s_shock_2;
else
    solverOutputs.shockPosition = [];
    solverOutputs.s_shock_1 = [];
    solverOutputs.s_shock_2 = [];
end
solverOutputs.s_injection = s_injection;
solverOutputs.s_chamber = s_chamber;
solverOutputs.s_throat = s_throat;
solverOutputs.s_exit = s_exit;
solverOutputs.optimalRatio = [];
solverOutputs.pressure_ideal = pressure_ideal;
solverOutputs.pressure_shockAtExit = pressure_shockAtExit;
solverOutputs.pressure_shockAtThroat = pressure_shockAtThroat;
solverOutputs.pressure_separation = pressure_separation;
solverOutputs.hasCondensation = gas.hasCondensation;
solverOutputs.expansionRatio_shockAtThroat = [];
solverOutputs.expansionRatio_optimal = [];
solverOutputs.expansionRatio_shockAtExit_1 = [];
solverOutputs.expansionRatio_shockAtExit_2 = [];
solverOutputs.expansionRatio_condensation = [];
solverOutputs.expansionRatio_separation = [];
solverOutputs.flowState = flowState;
end