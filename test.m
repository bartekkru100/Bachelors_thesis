import State.*
import Unit.*

% Test file to test parts of code

clc, clear, cleanup

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