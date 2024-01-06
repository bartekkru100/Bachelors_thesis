import Gas.*
import Unit.*
import PhaseDiagram.*
import State.*

% Test file to test parts of code


clc, clear, cleanup

atmo = Gas('air.yaml');
s_atmo = State(atmo);
s_atmo.pressure = oneatm * 0.1;
%s_atmo.pressure = 1.167594328952846e+07 * 1.2418;
setstate(atmo, s_atmo);
s_atmo = State(atmo);
expansionRatio = 460;
contractionRatio = 15;
Tolerance = 1e-12;
heatPower = 0;

gas = Gas('R1highT.yaml');
%gas = Gas();

%setstate(gas, 'Y', 'CH4:1,O2:3.6', 'T', 150, 'P', 300e5);
setstate(gas, 'Y', 'H2:1,O2:6.03', 'T', 150, 'P', 20.64e6);
%setstate(gas, 'Y', 'H2O:1', 'T', 1200, 'P', 10e5);
%setstate(gas, 'Y', 'POSF7688:1,O2: 2.36', 'T', 150, 'P', 26.7e6);
%setstate(gas, 'Y', 'POSF7688:1,O2: 2.63', 'T', 200, 'P', 24.52e6);
s_injection = State(gas);

setchamberconditions(gas, s_injection)
s_chamber = State(gas);

mass_rate_area_old = 0;
for i = 1 : 20
    setthroatconditions(gas);
    s_throat = State(gas);

    error_M = (mass_rate_area_old - s_throat.massFlowFlux) / mass_rate_area_old;

    if abs(error_M) < Tolerance
        break;
    end

    setinjectionconditions(gas, s_injection, contractionRatio);
    s_injection = State(gas);

    setchamberconditions(gas, s_injection, heatPower);
    s_chamber = State(gas);

    mass_rate_area_old = s_throat.massFlowFlux;
end

[maxExpansion_shock, s_shockAtExit] = checkfornormalshock(gas, s_throat, s_atmo);
[maxExpansion_separation, s_separation] = checkforseparation(gas, s_throat, s_atmo, 0);
setsupersonicexitconditions(gas, expansionRatio);
s_supersonicExit = State(gas);
[shockPosition, flowState, s_shock_1, s_shock_2] = shockposition(gas, s_atmo, s_supersonicExit, s_throat);
s_shockExit = State(gas);
s_shockStagnation = gas.stagnation;
setsubsonicexitconditions(gas, s_chamber, s_atmo);

specificImpulse = specificimpulse(gas, s_atmo);







disp("Injection: "+ "Velocity: " + s_injection.velocity         + " Entropy: " + s_injection.entropy        + " Total Energy: " + s_injection.totEnergy         + "Mass Flow Rate: " + s_injection.massFlowFlux * contractionRatio)
disp("Chamber:   "+ "Velocity: " + s_chamber.velocity           + " Entropy: " + s_chamber.entropy          + " Total Energy: " + s_chamber.totEnergy           + "Mass Flow Rate: " + s_chamber.massFlowFlux * contractionRatio)
disp("Throat:    "+ "Velocity: " + s_throat.velocity            + " Entropy: " + s_throat.entropy           + " Total Energy: " + s_throat.totEnergy            + "Mass Flow Rate: " + s_throat.massFlowFlux)
disp("Exit:      "+ "Velocity: " + s_supersonicExit.velocity    + " Entropy: " + s_supersonicExit.entropy   + " Total Energy: " + s_supersonicExit.totEnergy    + "Mass Flow Rate: " + s_supersonicExit.massFlowFlux * expansionRatio)
disp("ShockExit: "+ "Velocity: " + s_shockExit.velocity         + " Entropy: " + s_shockExit.entropy        + " Total Energy: " + s_shockExit.totEnergy         + "Mass Flow Rate: " + s_shockExit.massFlowFlux * shockPosition)
disp("ShockStag: "+ "Velocity: " + s_shockStagnation.velocity   + " Entropy: " + s_shockStagnation.entropy  + " Total Energy: " + s_shockStagnation.totEnergy   + "Mass Flow Rate: " + s_shockStagnation.massFlowFlux)
disp("Separat:   "+ "Velocity: " + s_separation.velocity        + " Entropy: " + s_separation.entropy       + " Total Energy: " + s_separation.totEnergy        + "Mass Flow Rate: " + s_separation.massFlowFlux * maxExpansion_separation)
disp("Gas:       "+ "Velocity: " + gas.velocity                 + " Entropy: " + gas.entropy                + " Total Energy: " + gas.totEnergy                 + "Mass Flow Rate: " + gas.massFlowFlux * expansionRatio)
%disp("Mach: " + s_shock_1.Mach + " p_sep/p_a: " + s_shock_1.pressure / s_atmo.pressure);

%{
phase = PhaseDiagram('D:\Uni\Inzynierka\Program\PhaseCurves\H2O')
point.temperature = 273.1;
stateofmatter().solid;
point.pressure = oneatm;
checkstateofmatter(phase, point)

unitArray = cell(1, 8);
mainFolder = dir("Units");
nSubFolders = size(mainFolder, 1);
for i = 3 : nSubFolders
    pathFolder = append(mainFolder(i).folder, '\' , mainFolder(i).name);
    subFolder = dir(pathFolder);
    nFiles = size(subFolder, 1);
    if nFiles > 2
        for j = 3 : nFiles
            if contains(subFolder(j).name, '.unit')
                pathFile = append(subFolder(j).folder, '\' , subFolder(j).name);
                unit = Unit(pathFile);
                unitArray{unit.type} = [unitArray{unit.type}, unit];
            end
        end
    end
end
%}
