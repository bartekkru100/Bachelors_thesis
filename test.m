import Gas.*
import Unit.*
import PhaseDiagram.*
import State.*

% Test file to test parts of code


clc, clear, cleanup

atmo = Gas('air.yaml');
s_atmo = State(atmo);
s_atmo.pressure = oneatm;
setstate(atmo, s_atmo);
s_atmo = State(atmo);
expansionRatio = 77.5;
contractionRatio = 10;
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

massFlowFlux_old = 0;
for i = 1 : 20

    [s_throat, s_stagnation] = gas.sonic;

    error_M = (massFlowFlux_old - s_throat.massFlowFlux) / massFlowFlux_old;
    if abs(error_M) < Tolerance
        break;
    end

    setinjectionconditions(gas, s_injection, contractionRatio);
    s_injection = State(gas);

    setchamberconditions(gas, s_injection, heatPower);
    s_chamber = State(gas);

    massFlowFlux_old = s_throat.massFlowFlux;
end

s_maxVelocity = gas.maxvelocity;

displaygasinfo(s_injection, "Injection", contractionRatio);
displaygasinfo(s_chamber, "Chamber", contractionRatio);
displaygasinfo(s_throat, "Throat", 1);



if s_stagnation.pressure > s_atmo.pressure
    [expansionRatio_shockAtExit_1, s_shockAtExit_1, s_shockAtExit_2, exitShockfound] = checkforshockatexit(gas, s_throat, s_atmo, s_maxVelocity, 'min');
    expansionRatio_shockAtThroat = checkforidealexpansion(gas, s_throat, s_atmo);

    if s_throat.pressure < s_atmo.pressure

        if expansionRatio < expansionRatio_shockAtThroat
            [s_exit, s_chamber] = setsubsonicexitconditions(gas, s_chamber, s_throat, s_atmo, expansionRatio, contractionRatio);
            flowState = "Fully subsonic";

        elseif expansionRatio <= expansionRatio_shockAtExit_1
            [expansionRatio_shockAtExit_2, s_shockAtExit_1, s_shockAtExit_2, exitShockfound] = checkforshockatexit(gas, s_throat, s_atmo, s_maxVelocity, max);

            if expansionRatio >= expansionRatio_shockAtExit_2
                expansionRatio = setsupersonicexitconditions(gas, s_throat, expansionRatio);
                flowState = "Fully supersonic";

            else
                [shockPosition, s_shock_1, s_shock_2] = shockposition(gas, s_atmo, s_supersonicExit, s_throat);
                flowState = "Normal shock inside the nozzle";

            end
        end
    else
        if expansionRatio <= expansionRatio_shockAtExit_1
            expansionRatio = setsupersonicexitconditions(gas, s_throat, expansionRatio);
            flowState = "Fully supersonic";

        else
            [shockPosition, s_shock_1, s_shock_2] = shockposition(gas, s_atmo, s_supersonicExit, s_throat);
            flowState = "Normal shock inside the nozzle";

        end
    end
else
    flowState = "No flow possible";
end
s_exit = State(gas);

if flowState == "Fully supersonic"
    if s_exit.pressure < s_atmo.pressure
    flowState = flowState + ", overexpanded";
    else
        flowState = flowState + ", underexpanded";
    end
end

displaygasinfo(s_exit, "exit", expansionRatio);

%setsubsonicexitconditions(gas, s_chamber, s_atmo);
specificImpulse = specificimpulse(gas, s_atmo);



function displaygasinfo(gas, name, areaRatio)
disp(name + ':' + ...
    " Velocity: " + gas.velocity + ...
    " Entropy: " + gas.entropy + ...
    " Total Energy: " + gas.totEnergy + ...
    " Mass Flow Rate: " + gas.massFlowFlux * areaRatio)
end

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
