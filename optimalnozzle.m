function solverOutputs = optimalnozzle(solverInputs)

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
setchamberconditions(gas, heatMass);
s_chamber = State(gas);
s_throat = gas.sonic;

s_stagnation = gas.stagnation;

expansionRatio_shockAtThroat = [];
expansionRatio_shockAtExit_2 = [];
expansionRatio_separation = [];
setstate(gas, s_throat);
expansionRatio_condensation = findcondensation(gas);
setstate(gas, s_throat);
expansionRatio_optimal = findidealexpansion(gas, atmo);


[expansionRatio_shockAtExit_1, shockFound] = findshockatexit(gas, atmo, 'min');
if atmo.pressure > s_stagnation.pressure
    flowState = "Stagation pressure below ambient, flow impossible";

    s_injection = State(setstate(gas, 'P', atmo.pressure, 'T', atmo.temperature, 'velocity', 0));
    s_chamber = s_injection;
    s_throat = s_injection;
    s_exit = s_injection;
elseif atmo.pressure > s_throat.pressure
    if shockFound
        flowState = "Throat pressure below ambient, optimal nozzle subsonic";
        setstate(gas, s_throat);
        expansionRatio_shockAtExit_2 = findshockatexit(gas, atmo, 'max');
        expansionRatio_separation = 1;
        expansionRatio_shockAtThroat = expansionRatio_optimal;

        [s_exit, s_throat, s_chamber] = setsubsonicexitconditions(gas, atmo, expansionRatio_optimal, contractionRatio);
    else
        flowState = "Throat pressure below ambient, optimal nozzle subsonic, fully supersonic flow impossible";

        expansionRatio_shockAtThroat = expansionRatio_optimal;
        expansionRatio_shockAtExit_1 = [];

        setstate(gas, s_throat);
        [s_exit, s_throat, s_chamber] = setsubsonicexitconditions(gas, atmo, expansionRatio_optimal, contractionRatio);
    end
else
    flowState = "Optimal nozzle supersonic";

    expansionRatio_shockAtExit_2 = [];
    setstate(gas, s_throat);
    expansionRatio_separation = findseparation(gas, atmo, separationTolerance);
    setpressureisentropic(gas, atmo.pressure, s_throat);
    s_exit = State(gas);
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

solverOutputs.A_chamber = [];
solverOutputs.A_throat = [];
solverOutputs.A_exit = [];
solverOutputs.D_chamber = [];
solverOutputs.D_throat = [];
solverOutputs.D_exit = [];
solverOutputs.massFlow = [];
solverOutputs.heatPower = [];

solverOutputs.shockPosition = [];
solverOutputs.s_shock_1 = [];
solverOutputs.s_shock_2 = [];

solverOutputs.s_injection = s_injection;
solverOutputs.s_chamber = s_chamber;
solverOutputs.s_throat = s_throat;
solverOutputs.s_exit = s_exit;
solverOutputs.optimalRatio = [];
solverOutputs.pressure_ideal = [];
solverOutputs.pressure_shockAtExit = [];
solverOutputs.pressure_shockAtThroat = [];
solverOutputs.pressure_separation = [];
solverOutputs.hasCondensation = gas.hasCondensation;
solverOutputs.expansionRatio_shockAtThroat = expansionRatio_shockAtThroat;
solverOutputs.expansionRatio_optimal = expansionRatio_optimal;
solverOutputs.expansionRatio_shockAtExit_1 = expansionRatio_shockAtExit_1;
solverOutputs.expansionRatio_shockAtExit_2 = expansionRatio_shockAtExit_2;
solverOutputs.expansionRatio_condensation = expansionRatio_condensation;
solverOutputs.expansionRatio_separation = expansionRatio_separation;
solverOutputs.flowState = flowState;
end