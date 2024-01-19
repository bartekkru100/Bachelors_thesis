classdef PhaseDiagram

    % This class stores data about a substance's phase diagram

    properties(Access = private)
        line_lg % liquid-gas curve
        line_sg % solid-gas curve
        line_ls % liquid-solid curve
        src; % source directory
        species; % name of the species
    end

    methods(Access = public)
        
%==========================================================================

        % The construction takes 3 curve files from a chosen directory

        function phaseDiagram = PhaseDiagram(source)
            import PhaseDiagram.*

            phaseDiagram.src = source;
            sourceFolder = dir(source);
            test = (strfind(source, '\'));
            phaseDiagram.species = source(test(end) + 1 : end);
            nFiles = size(sourceFolder, 1);
            for j = 3 : nFiles
                fileName = lower(sourceFolder(j).name);
                if ismember(fileName, ['lg.txt', 'sg.txt', 'ls.txt'])
                    Table = readtable(append(sourceFolder(j).folder, '\', sourceFolder(j).name));
                    if fileName == 'lg.txt'
                        phaseDiagram.line_lg = table2line(Table);
                    end
                    if fileName == 'sg.txt'
                        phaseDiagram.line_sg = table2line(Table);
                    end
                    if fileName == 'ls.txt'
                        phaseDiagram.line_ls = table2line(Table);
                    end
                end
            end
        end

%==========================================================================

        % Getters

        function src = source(phaseDiagram)
            src = phaseDiagram.src;
        end

        function name = speciesName(phaseDiagram)
            name = phaseDiagram.species;
        end

%==========================================================================

        % This checks the what a state of matter is at a given point on P-T
        % diagram

        function stateOfMatter = checkstateofmatter(phaseDiagram, point)
            import PhaseDiagram.*

            isCondensed = false;
            if isbounded(point.pressure, phaseDiagram.line_lg.pressure)
                if point.temperature < getTemperature(point.pressure, phaseDiagram.line_lg)
                    isCondensed = true;
                end
            elseif isbounded_max(point.pressure, phaseDiagram.line_sg.pressure)
                if point.temperature < getTemperature(point.pressure, phaseDiagram.line_sg)
                    isCondensed = true;
                end
            end
            if isCondensed
                if point.temperature > getTemperature(point.pressure, phaseDiagram.line_ls)
                    stateOfMatter = stateofmatter().liquid;
                else
                    stateOfMatter = stateofmatter().solid;
                end
            else
                stateOfMatter = stateofmatter().gas;
            end
        end
    end

%==========================================================================

    % Enum

    methods (Static, Access = public)
        function stateOfMatter = stateofmatter()
            stateOfMatter = struct('gas', 1, 'liquid', 2, 'solid', 3);
        end

        function stateOfMatter = tostateofmatter(str)
            switch lower(convertStringsToChars(lower(str)))
                case 'gas'
                    stateOfMatter = 1;
                case 'liquid'
                    stateOfMatter = 2;
                case 'solid'
                    stateOfMatter = 3;
            end
        end
    end

    methods (Static, Access = private)
        
%==========================================================================

        % Some utility functions for inside use

        function line = table2line(Table)
            line.temperature = table2array(Table(:,1));
            line.pressure = table2array(Table(:,2));
        end

        function isBounded = isbounded(value, array)
            if value >= min(array) && value <= max(array)
                isBounded = true;
            else
                isBounded = false;
            end
        end

        function isBounded = isbounded_min(value, array)
            if value >= min(array)
                isBounded = true;
            else
                isBounded = false;
            end
        end

        function isBounded = isbounded_max(value, array)
            if value <= max(array)
                isBounded = true;
            else
                isBounded = false;
            end
        end

        function pressure = getPressure(temperature, line)
            pressure = interp1(line.temperature, line.pressure, temperature, "makima");
        end

        function temperature = getTemperature(pressure, line)
            temperature = interp1(line.pressure, line.temperature, pressure, "makima");
        end
    end
end