function solverOutputs = engineanalysis(solverInputs)

import Gas.*
import State.*

propellantArray = solverInputs.propellantArray;
gas = solverInputs.gas;
ratios = solverInputs.ratios;
atmo = solverInputs.atmo;
expansionRatio = solverInputs.expansionRatio;
contractionRatio = solverInputs.contractionRatio;
heatMass = solverInputs.heatMass;
thrust = solverInputs.thrust;
separationTolerance = solverInputs.separationTolerance;
mass_mole = solverInputs.mass_mole;

gas = combine(gas, propellantArray, ratios, mass_mole);
s_injection = State(gas);
setchamberconditions(gas, heatMass);
setarearatioisentropic(gas, contractionRatio, "sub");
s_chamber = State(gas);
s_throat = gas.sonic;
s_exit = setsupersonicexitconditions(gas, expansionRatio);
s_stagnation = gas.stagnation;

pressure_shockAtThroat = findshockatthroat_pressure(gas, expansionRatio);
pressure_shockAtExit = findshockatexit_pressure(gas, s_exit, expansionRatio);
pressure_separation  = findseparation_pressure(gas, s_exit, separationTolerance);
pressure_ideal = findidealexpansion_pressure(gas, s_exit);
%pressure_condensation = findcondensation_pressure(gas);

shockFound = false;
separationFound = false;

flowState = Enum("No flow possible", "Fully subsonic", "Normal shock inside the nozzle", "Fully supersonic, underexpanded", "Fully supersonic, overexpanded");

if atmo.pressure > s_stagnation.pressure
    flowState.current = flowState.value.no_flow_possible;

    s_injection = State(atmo);
    s_chamber = s_injection;
    s_throat = s_injection;
    s_exit = s_injection;
else
    if atmo.pressure > pressure_shockAtThroat
        [s_exit, s_throat, s_chamber] = setsubsonicexitconditions(gas, atmo, expansionRatio, contractionRatio);
        flowState.current = flowState.value.fully_subsonic;
    else
        if atmo.pressure > pressure_shockAtExit
            [shockPosition, s_shock_1, s_shock_2, s_exit] = shockposition(gas, atmo, s_exit);
            flowState.current = flowState.value.normal_shock_inside_the_nozzle;
            shockFound = true;
        else
            if atmo.pressure >= pressure_ideal
                flowState.current = flowState.value.fully_supersonic_overexpanded;
                if atmo.pressure > pressure_separation
                    separationFound = true;
                end
            else
                flowState.current = flowState.value.fully_supersonic_underexpanded;
            end
        end
    end
end

effectiveExhaustVelocity = effectiveexhaustvelocity(gas, atmo);
massFlow = thrust / effectiveExhaustVelocity;
A_chamber = massFlow / s_chamber.massFlowFlux;
A_throat = massFlow / s_throat.massFlowFlux;
A_exit = massFlow / s_exit.massFlowFlux;
heatPower = heatMass * massFlow;
D_chamber = sqrt(4 * A_chamber / pi);
D_throat = sqrt(4 * A_throat / pi);
D_exit = sqrt(4 * A_exit / pi);

solverOutputs.s_injection = s_injection;
solverOutputs.s_chamber = s_chamber;
solverOutputs.s_throat = s_throat;
solverOutputs.s_exit = s_exit;
solverOutputs.pressure_ideal = pressure_ideal;
solverOutputs.pressure_shockAtExit = pressure_shockAtExit;
solverOutputs.pressure_shockAtThroat = pressure_shockAtThroat;
solverOutputs.pressure_separation = pressure_separation;
solverOutputs.hasCondensation = gas.hasCondensation;
solverOutputs.effectiveExhaustVelocity = effectiveExhaustVelocity;
solverOutputs.A_chamber = A_chamber;
solverOutputs.A_throat = A_throat;
solverOutputs.A_exit = A_exit;
solverOutputs.D_chamber = D_chamber;
solverOutputs.D_throat = D_throat;
solverOutputs.D_exit = D_exit;
solverOutputs.massFlow = massFlow;
solverOutputs.heatPower = heatPower;
solverOutputs.flowState = flowState;
solverOutputs.shockFound = shockFound;
solverOutputs.separationFound = separationFound;
if shockFound
    solverOutputs.shockPosition = shockPosition;
    solverOutputs.s_shock_1 = s_shock_1;
    solverOutputs.s_shock_2 = s_shock_2;
else
    solverOutputs.shockPosition = 1;
    solverOutputs.s_shock_1 = s_exit;
    solverOutputs.s_shock_2 = s_exit;
end
end