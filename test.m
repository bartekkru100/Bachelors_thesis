import Gas.*
addpathall(pwd, ".git", true);
%adddir('D:\Uni\Inzynierka\Program')
clc, clear, cleanup

prop = Gas('GRI30')
setsupersonicexitconditions(prop, prop.sonic, 10)


atmo = Gas('air');
setstate(atmo, 'P', 0.00001, 'T', atmo.temperature);

propellantArray(2) = Gas('GRI30');
setstate(propellantArray(2), 'Y', 'O2:1', 'T', 150, 'P', 20.64e6);
propellantArray(1) = Gas('GRI30');
setstate(propellantArray(1), 'Y', 'H2:1', 'T', 150, 'P', 20.64e6);
ratios = [1, 6.03];
gas = combine('new', propellantArray, ratios, 'mass');

contractionRatio = 15;
expansionRatio = 77.5;
thrust = 1.86e6;

temperature_chamber_max = 4000;
temperature_throat_max = 3000;
heatPowerArea = 0;

solverInputs.gas = gas;
solverInputs.atmo = atmo;
solverInputs.propellantArray = propellantArray;
solverInputs.ratios = ratios;
solverInputs.contractionRatio = contractionRatio;
solverInputs.expansionRatio = expansionRatio;
solverInputs.thrust = thrust;
solverInputs.heatPowerArea = heatPowerArea;
solverInputs.temperature_chamber_max = temperature_chamber_max;
solverInputs.temperature_throat_max = temperature_throat_max;

engineanalysis(solverInputs);
%optimalmixture(solverInputs);