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
            import Unit.*
            if nargin == 1
                unit = unit.readfromfile(name);
            else
                unit.name = name;
                unit.type = unitTypes().fromstr(type);
                unit.symbol = symbol;
                unit.multiplier = multiplier;
                if nargin == 5
                    unit.offset = offset;
                else
                    unit.offset = 0;
                end
            end

        end

        function unitValue = apply(unit, value, dimension)
            arguments
                unit;
                value;
                dimension = 1;
            end
            unitValue = (value + unit.offset) / unit.multiplier ^ dimension;
        end

        function type = gettype(this)
            type = this.type;
        end

        function name = getname(this)
            name = this.name;
        end

        function symbol = getsymbol(this, dimension)
            if nargin == 1 || dimension == 1
                symbol = convertStringsToChars(this.symbol);
            else
                symbol = convertStringsToChars(this.symbol) + "^" + num2str(dimension);
            end
        end

        function multiplier = getmultiplier(this)
            multiplier = this.multiplier;
        end

        function offset = getoffset(this)
            offset = this.offset;
        end

        function value = remove(unit, unitValue, dimension)
            arguments
                unit;
                unitValue;
                dimension = 1;
            end
            value = unitValue * unit.multiplier ^ (1 / dimension) - unit.offset;
        end

%==========================================================================

        % This loads unit info from a file

        function unit = readfromfile(unit, fileName)
            import Unit.*
            
            file = fopen(fileName, 'r');
            unit.type = unitTypes().fromstr(fscanf(file, 'Type: %s\n'));
            unit.name = fscanf(file, 'Name: %s\n');
            unit.symbol = fscanf(file, 'Symbol: %s\n');
            unit.multiplier = fscanf(file, 'Multiplier: %f\n');
            unit.offset = fscanf(file, 'Offset: %f');
            fclose(file);
        end

    end

%==========================================================================

    % Enum of unit types

    methods (Static, Access = public)
        
        function unitType = unitTypes()
            import Enum.*
            unitType = Enum('length', 'mass', 'speed', 'force', 'pressure', 'temperature', 'energy', 'power', 'ratio');
        end

    end

end
