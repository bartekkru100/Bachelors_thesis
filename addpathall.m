function addpathall(currentFolder, exclusions, isRecursive)
if nargin < 3
    isRecursive = false;
end
subFolders = ls(currentFolder);
nSubFolders = size(subFolders, 1);
for i = 3 : nSubFolders
    subFolder = [currentFolder, '\', strtrim(subFolders(i,:))];
    if nargin < 2 || ~ismember(convertCharsToStrings(strtrim(subFolders(i,:))), exclusions)
        if isfolder(subFolder)
            addpath(subFolder);
            if isRecursive
                addpathall(subFolder, string.empty, isRecursive);
            end
        end
    end
end
end