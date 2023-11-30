classdef Unit

    properties (Access = public)

        type;
        name;
        symbol;
        multiplier;
        offset;

    end

    methods (Access = public)

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

        function unit = readfromfile(unit, fileName)
            file = fopen(fileName, 'r');
            unit.type = fscanf(file, 'Type: %s\n');
            unit.name = fscanf(file, 'Name: %s\n');
            unit.symbol = fscanf(file, 'Symbol: %s\n');
            unit.multiplier = fscanf(file, 'Multiplier: %f\n');
            unit.offset = fscanf(file, 'Offset: %f');
            fclose(file);
        end

    end

    methods (Static)

        function aunitValue = applyUnit(value, unit)
            import this.*
            aunitValue = unit.apply(value);
        end

        function value = removeUnit(unitValue, unit)
            import this.*
            value = unit.remove(unitValue);
        end

    end

end
