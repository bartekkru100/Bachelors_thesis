import Gas.*
import Unit.*
import State.*

% Test file to test parts of code


clc, clear, cleanup
soot = Solution('graphite.yaml');
gas = GRI30;
exhaust = Mixture({gas, 1.0; soot, 0.0});
set(gas, 'P', 10 * oneatm)

%{
atmo = Gas('air.yaml');
s_atmo = State(atmo);
s_atmo.pressure = oneatm;
setstate(atmo, s_atmo);
s_atmo = State(atmo);
expansionRatio = 100000;
contractionRatio = 10;

gas = Gas();

setstate(gas, 'Y', 'CH4:1,O2:3.6', 'T', 150, 'P', 300e5);
s_injection = State(gas);
 
setchamberconditions(gas, s_injection)
s_chamber = State(gas);

setthroatconditions(gas);
s_throat = State(gas);

setsupersonicexitconditions(gas, expansionRatio);
s_exit = State(gas);

[shock, flowState] = shockposition(gas, s_atmo, s_exit, s_throat);
s_exit_s = State(gas);
disp("Velocity: " + s_injection.velocity + " Entropy:" + s_injection.entropy + " Total Energy: " + s_injection.totEnergy + "Mass Flow Rate: " + s_injection.massFlowFlux * contractionRatio)
disp("Velocity: " + s_chamber.velocity + " Entropy: " + s_chamber.entropy + " Total Energy: " + s_chamber.totEnergy + "Mass Flow Rate: " + s_chamber.massFlowFlux * contractionRatio)
disp("Velocity: " + s_throat.velocity + " Entropy: " + s_throat.entropy + " Total Energy: " + s_throat.totEnergy + "Mass Flow Rate: " + s_throat.massFlowFlux)
disp("Velocity: " + s_exit.velocity + " Entropy: " + s_exit.entropy + " Total Energy: " + s_exit.totEnergy + "Mass Flow Rate: " + s_exit.massFlowFlux * expansionRatio)
disp("Velocity: " + gas.velocity + " Entropy: " + gas.entropy + " Total Energy: " + gas.totEnergy + "Mass Flow Rate: " + gas.massFlowFlux * shock)
%}


%{

%unitTypes = {'length', 'mass', 'speed', 'force', 'pressure', 'temperature', 'energy', 'power'};
unitArray = dictionary();
unitArray("length") = {[Unit("Units\length\meter.unit"), Unit("Units\length\milimeter.unit")]}
unitArray("length") = {cell2mat(unitArray("length")) + Unit("Units\length\meter.unit")}

nUnits = dictionary();
%nUnits('length') = 1
mainFolder = dir("Units");
nSubFolders = size(mainFolder, 1);
for i = 3 : nSubFolders;
    pathFolder = append(mainFolder(i).folder, '\' , mainFolder(i).name);
    subFolder = dir(pathFolder);
    nFiles = size(subFolder, 1);
    if nFiles > 2
        for j = 1 : nFiles
            if contains(subFolder(j).name, '.unit')
                pathFile = append(subFolder(j).folder, '\' , subFolder(j).name);
                unitType = Unit(pathFile).type
                if isConfigured(nUnits)
                    nUnits(unitType) = nUnits(unitType) + 1;
                else
                    nUnits(unitType) = 1;
                end
                unitArray(unitType) = {[Unit(pathFile)]};
            end
        end
    end
end
%}