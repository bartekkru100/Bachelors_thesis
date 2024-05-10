addpathall(pwd, ".git", true);
addpath('D:\Uni\Inzynierka\Program')
clc, clear%, cleanup
import Gas.*

atmo = Gas('air');
gas = Gas('GRI30');
i = 0;
tic
for heatMass = logspace(5,7,10)
    i = i + 1;
    setstate(gas, 'Y', 'O2:1', 'T', 50, 'P', 20e5, 'velocity', 0);
    propellantArray(1) = State(gas);
    setstate(gas, 'Y', 'H2:1', 'T', 50, 'P', 20e5, 'velocity', 0);
    propellantArray(2) = State(gas);
    setstate(atmo, 'P', oneatm, 'T', atmo.temperature);

    solverInputs.atmo = atmo;
    solverInputs.gas = gas;
    solverInputs.propellantArray = propellantArray;
    solverInputs.ratios = [1, 6];
    solverInputs.contractionRatio = 5;
    solverInputs.expansionRatio = 69;
    solverInputs.thrust = 1.0331e6;
    solverInputs.heatMass = heatMass;
    solverInputs.pressure_injection = 200e5;
    solverInputs.temperature_injection = 50;
    solverInputs.temperature_chamber_max = 2000;
    solverInputs.temperature_throat_max = 2000;
    solverInputs.separationTolerance = 0;
    solverInputs.mass_mole = "mass";
    solverInputs.searchedIndicator = "specificImpulse";
    % solverOutputs = engineanalysis(solverInputs);
    % solverOutputs = optimalnozzle(solverInputs);
    solverOutputs = optimalmixture(solverInputs);
    x(i) = heatMass;
    y(i) = solverOutputs.optimalRatio;
end
toc
semilogx(x, y)
