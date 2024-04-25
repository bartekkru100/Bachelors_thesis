addpathall(pwd, ".git", true);
adddir('D:\Uni\Inzynierka\Program')
clc, clear, cleanup

import Gas.*

atmo = Gas('air');
setstate(atmo, 'P', oneatm, 'T', atmo.temperature);
%sol = Solution('nasa_gas.yaml')
gas = Gas('gri30');
setstate(gas, 'Y', 'O2:1', 'T', 50, 'P', 1.5e5);
propellantArray(1) = State(gas);
setstate(gas, 'Y', 'CH4:1', 'T', 50, 'P', 1.5e5);
propellantArray(2) = State(gas);
ratios = [6.03, 1];
mass_mole = "mass";
contractionRatio = 5;
expansionRatio = 20;
thrust = 1.0331e6;

i= 0;
for pressure = 1e5 : 1e3 : 2e5
    i = i + 1;
setstate(gas, 'Y', 'H2:1', 'T', 50, 'P', 2e5, 'velocity', 0);
propellantArray(1) = State(gas);
setstate(gas, 'Y', 'O2:1', 'T', 50, 'P', 2e5, 'velocity', 0);
propellantArray(2) = State(gas);
setstate(atmo, 'P', pressure, 'T', atmo.temperature);

temperature_chamber_max = 4000;
temperature_throat_max = 3000;
separationTolerance = 0;
heatMass = 0;
solverInputs.atmo = atmo;
solverInputs.gas = gas;
solverInputs.propellantArray = propellantArray;
solverInputs.ratios = ratios;
solverInputs.contractionRatio = contractionRatio;
solverInputs.expansionRatio = expansionRatio;
solverInputs.thrust = thrust;
solverInputs.heatMass = heatMass;
solverInputs.temperature_chamber_max = temperature_chamber_max;
solverInputs.temperature_throat_max = temperature_throat_max;
solverInputs.separationTolerance = separationTolerance;
solverInputs.mass_mole = mass_mole;
solverOutputs = engineanalysis(solverInputs);
%solverOutputs = optimalmixture(solverInputs);
% solverOutputs = optimalnozzle(solverInputs);
% log(i, 1) = expansionRatio;
% log(i, 2) = solverOutputs.pressure_ideal;
% log(i, 3) = solverOutputs.pressure_separation;
% log(i, 4) = solverOutputs.pressure_shockAtExit;
% log(i, 5) = solverOutputs.pressure_shockAtThroat;
end
% plot(log(:, 1),log(:, 2));
% hold on;
% plot(log(:, 1),log(:, 3));
% hold on;
% plot(log(:, 1),log(:, 4));
% hold on;
% plot(log(:, 1),log(:, 5));