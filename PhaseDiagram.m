classdef PhaseDiagram
    properties(Access = private)
        line_lg
        line_sg
        line_ls
        src;
        species;
    end

    methods(Access = public)
        function phaseDiagram = PhaseDiagram(source)
            import PhaseDiagram.*

            sourceFolder = dir(source);
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

        function src = source(phaseDiagram)
            src = phaseDiagram.src
        end

        function name = speciesName(phaseDiagram)
            name = phaseDiagram.species
        end

        function stateOfMatter = checkstateofmatter(phaseDiagram, point)
            import PhaseDiagram.*

            isCondensed = false;
            if isbounded(point.temperature, phaseDiagram.line_lg.temperature)
                if point.pressure > getPressure(point.temperature, phaseDiagram.line_lg)
                    getPressure(point.temperature, phaseDiagram.line_lg)
                    isCondensed = true;
                end
            elseif isbounded(point.temperature, phaseDiagram.line_sg.temperature)
                if point.pressure > getPressure(point.temperature, phaseDiagram.line_sg)
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

        function pressure = getPressure(temperature, line)
            pressure = interp1(line.temperature, line.pressure, temperature, "spline");
        end

        function temperature = getTemperature(pressure, line)
            temperature = interp1(line.pressure, line.temperature, pressure, "spline");
        end
    end
end