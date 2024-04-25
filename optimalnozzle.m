function solverOutputs = optimalnozzle(solverInputs)

import Gas.*
import State.*

propellantArray = solverInputs.propellantArray;
gas = solverInputs.gas;
ratios = solverInputs.ratios;
atmo = solverInputs.atmo;
expansionRatio = solverInputs.expansionRatio;
contractionRatio = solverInputs.contractionRatio;
heatMass = solverInputs.heatMass;
separationTolerance = solverInputs.separationTolerance;
mass_mole = solverInputs.mass_mole;

gas = combine(gas, propellantArray, ratios, mass_mole);

s_injection = State(gas);
setchamberconditions(gas, heatMass);
%setarearatioisentropic(gas, contractionRatio, "sub");
s_chamber = State(gas);
s_throat = gas.sonic;

s_stagnation = gas.stagnation;


expansionRatio_condensation = findcondensation(gas);

expansionRatio_shockAtThroat = findidealexpansion(gas, atmo);
[expansionRatio_shockAtExit_1, shockFound] = findshockatexit(gas, atmo, 'min');
if shockFound
    if atmo.pressure > s_throat.pressure
        expansionRatio_shockAtExit_2 = findshockatexit(gas, atmo, 'max');
        expansionRatio_separation = 1;
    else
        expansionRatio_shockAtExit_2 = 1;
        expansionRatio_separation = findseparation(gas, atmo, separationTolerance);
    end
else
        expansionRatio_shockAtExit_1 = 1;
        expansionRatio_shockAtExit_2 = 1;
        expansionRatio_separation = 1;
end




s_exit = setsupersonicexitconditions(gas, expansionRatio);

solverOutputs.s_injection = s_injection;
solverOutputs.s_chamber = s_chamber;
solverOutputs.s_throat = s_throat;
solverOutputs.s_exit = s_exit;
solverOutputs.expansionRatio_shockAtThroat = expansionRatio_shockAtThroat;
solverOutputs.expansionRatio_shockAtExit_1 = expansionRatio_shockAtExit_1;
solverOutputs.expansionRatio_shockAtExit_2 = expansionRatio_shockAtExit_2;
solverOutputs.expansionRatio_condensation = expansionRatio_condensation;
solverOutputs.expansionRatio_separation = expansionRatio_separation;
end