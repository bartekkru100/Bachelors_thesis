classdef Unit

    % This class stores unit data for conversions

    properties (Access = public)

        type;
        name;
        symbol;
        multiplier;
        offset; % Used for temperatures

    end

    methods (Access = public)

%==========================================================================

        % The constructor

        function unit = Unit(name, type, symbol, multiplier, offset)
            if nargin == 1
                unit = unit.readfromfile(name);
            else
                unit.name = name;
                unit.type = type;
                unit.symbol = symbol;
                unit.multiplier = multiplier;
                if nargin == 5
                    unit.offset = offset;
                else
                    unit.offset = 0;
                end
            end

        end

        function unitValue = apply(unit, value)
            unitValue = (value + unit.offset) / unit.multiplier;
        end

        function Value = remove(unit, unitValue)
            Value = unitValue * unit.multiplier - unit.offset;
        end

%==========================================================================

        % This loads unit info from a file

        function unit = readfromfile(unit, fileName)
            import Unit.*
            
            file = fopen(fileName, 'r');
            unit.type = toUnitType(fscanf(file, 'Type: %s\n'));
            unit.name = fscanf(file, 'Name: %s\n');
            unit.symbol = fscanf(file, 'Symbol: %s\n');
            unit.multiplier = fscanf(file, 'Multiplier: %f\n');
            unit.offset = fscanf(file, 'Offset: %f');
            fclose(file);
        end

    end

%==========================================================================

    % Enum of unit types

    methods (Static)
        
        function type = unitType()
        type = struct('length', 1, 'mass', 2, 'speed', 3, 'force', 4, 'pressure', 5, 'temperature', 6, 'energy', 7, 'power', 8);
        end

        function type = tounittype(str)
            switch lower(convertStringsToChars(lower(str)))
                case 'length'
                    type = 1;
                case 'mass'
                    type = 2;
                case 'speed'
                    type = 3;
                case 'force'
                    type = 4;
                case 'pressure'
                    type = 5;
                case 'temperature'
                    type = 6;
                case 'energy'
                    type = 7;
                case 'power'
                    type = 8;
            end
        end
        
        function unitValue = applyUnit(value, unit)
            import this.*
            unitValue = unit.apply(value);
        end

        function value = removeUnit(unitValue, unit)
            import this.*
            value = unit.remove(unitValue);
        end

    end

end
