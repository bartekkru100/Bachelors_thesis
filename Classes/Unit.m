classdef Unit

    % This class stores unit data for conversions

    properties (Access = public)

        type;
        name;
        symbol;
        multiplier;
        dimension;
        offset; % Used for temperatures

    end

    properties (Access = public)

        subUnits;
        symbol_full;
        latex_full;
        multiplier_full;
        is_compound;

    end

    methods (Access = public)

        %==========================================================================

        % The constructor

        function unit = Unit(name, type, symbol, multiplier, offset, dimension)
            import Unit.*
            if nargin == 1
                fileName = name;
                file = fopen(fileName, 'r');
                unit.type = unitTypes().fromstr(fscanf(file, 'Type: %s\n'));
                unit.name = convertCharsToStrings(fscanf(file, 'Name: %s\n'));
                unit.symbol = convertCharsToStrings(fscanf(file, 'Symbol: %s\n'));
                unit.multiplier = fscanf(file, 'Multiplier: %f\n');
                unit.offset = fscanf(file, 'Offset: %f');
                unit.dimension = 1;
                unit.is_compound = false;
                fclose(file);
            elseif nargin == 2
                subUnits = name;
                dimension = type;
                unit.subUnits = subUnits;
                unit.name = arrayfun(@(subUnit) subUnit.name, subUnits);
                unit.type = arrayfun(@(subUnit) subUnit.type, subUnits);
                unit.symbol = arrayfun(@(subUnit) subUnit.symbol, subUnits);
                unit.multiplier = arrayfun(@(subUnit) subUnit.multiplier, subUnits);
                if dimension == 1
                    unit.offset = subUnits.offset;
                else
                    unit.offset = 0;
                end
                unit.dimension = dimension;
                unit.is_compound = true;
            else
                if class(type) == "string" || class(type) == "char"
                    type = unitTypes().fromstr(type);
                end
                unit.name = convertCharsToStrings(name);
                unit.type = convertCharsToStrings(type);
                unit.symbol = convertCharsToStrings(symbol);
                unit.multiplier = convertCharsToStrings(multiplier);
                unit.offset = convertCharsToStrings(offset);
                unit.dimension = convertCharsToStrings(dimension);
                unit.is_compound = false;
            end
            unit.multiplier_full = prod(unit.multiplier .^ unit.dimension);

            nominatorSymbols = find(unit.dimension > 0);
            denominatorSymbols = find(unit.dimension < 0);


            unit.symbol_full = "";
            if isempty(nominatorSymbols)
                unit.symbol_full = "1";
            else
                firstIteration = true;
                for nominatorSymbol = nominatorSymbols
                    if firstIteration
                        firstIteration = false;
                    else
                        unit.symbol_full = unit.symbol_full + " ";
                    end
                    unit.symbol_full = unit.symbol_full + unit.symbol(nominatorSymbol);
                    if unit.dimension(nominatorSymbol) ~= 1
                        unit.symbol_full = unit.symbol_full + "^" + abs(unit.dimension(nominatorSymbol));
                    end
                end
            end
            if ~isempty(denominatorSymbols)
                unit.symbol_full = "(" + unit.symbol_full + ")/(";
                firstIteration = true;
                for denominatorSymbol = denominatorSymbols
                    if firstIteration
                        firstIteration = false;
                    else
                        unit.symbol_full = unit.symbol_full + " ";
                    end
                    unit.symbol_full = unit.symbol_full + unit.symbol(denominatorSymbol);
                    if unit.dimension(denominatorSymbol) ~= 1
                        unit.symbol_full = unit.symbol_full + "^" + abs(unit.dimension(denominatorSymbol));
                    end
                end
                unit.symbol_full = unit.symbol_full + ")";
            end

            unit.latex_full = "";
            nominatorLatex = "";
            denominatorLatex = "";
            if isempty(nominatorSymbols)
                nominatorLatex = "1";
            else
                firstIteration = true;
                for nominatorSymbol = nominatorSymbols
                    if firstIteration
                        firstIteration = false;
                    else
                        nominatorLatex = nominatorLatex + " ";
                    end
                    nominatorLatex = nominatorLatex + unit.symbol(nominatorSymbol);
                    if unit.dimension(nominatorSymbol) ~= 1
                        nominatorLatex = nominatorLatex + "^{" + unit.dimension(nominatorSymbol) + "}";
                    end
                end
            end
            if ~isempty(denominatorSymbols)
                firstIteration = true;
                for denominatorSymbol = denominatorSymbols
                    if firstIteration
                        firstIteration = false;
                    else
                        denominatorSymbol = denominatorSymbol + " ";
                    end
                    denominatorLatex = denominatorLatex + unit.symbol(denominatorSymbol);
                    if unit.dimension(denominatorSymbol) ~= -1
                        denominatorLatex = denominatorLatex + "^{" + abs(unit.dimension(denominatorSymbol)) + "}";
                    end
                end
                unit.latex_full = "\frac{" + nominatorLatex + "}{" + denominatorLatex + "}";
            else
                unit.latex_full = nominatorLatex;
            end

        end

        function unitValue = apply(unit, value, dimension)
            arguments
                unit;
                value;
                dimension = 1;
            end
            unitValue = (value + unit.offset) / unit.multiplier_full ^ (dimension);
        end

        function value = remove(unit, unitValue, dimension)
            arguments
                unit;
                unitValue;
                dimension = 1;
            end
            value = unitValue * unit.multiplier_full ^ (1 / (dimension)) - unit.offset;
        end

        function type = gettype(this)
            type = this.type;
        end

        function name = getname(this)
            name = this.name;
        end

        function symbol = getsymbol(this, dimension)
            arguments
                this;
                dimension = 1;
            end
            symbol = this.symbol_full;
            if dimension ~= 1
                if size(this.dimension, 2) > 1
                    symbol = "(" + symbol + ")";
                end
                symbol = symbol + "^" + dimension
            end
        end

        function symbol = getlatex(this, dimension)
            arguments
                this;
                dimension = 1;
            end
            symbol = this.latex_full;
            if dimension ~= 1
                if size(this.dimension, 2) > 1
                    if ~isempty(this.dimension(this.dimension < 0))
                        symbol = "\left(" + symbol + "\right)";
                    else
                        symbol = "(" + symbol + ")";
                    end
                end
                symbol = symbol + "^{" + dimension + "}";
            end
            symbol = strrep(symbol, "%", "\%");
            symbol = strrep(symbol, "‰", "\permil");
            symbol = "$" + symbol + "$";
        end

        function multiplier = getmultiplier(this)
            multiplier = this.multiplier;
        end

        function offset = getoffset(this)
            offset = this.offset;
        end

    end

    methods (Static, Access = public)

        %==========================================================================

        % Enum of unit types

        function unitType = unitTypes()
            import Enum.*
            unitType = Enum("length", "mass", "speed", "force", "pressure", "temperature", "energy", "power", "ratio", "time");
        end

        %==========================================================================

        % This loads unit info from a file

        function unit = readfromfile(fileName)
            import Unit.*

            file = fopen(fileName, 'r');
            type = fscanf(file, 'Type: %s\n');
            name = fscanf(file, 'Name: %s\n');
            symbol = fscanf(file, 'Symbol: %s\n');
            multiplier = fscanf(file, 'Multiplier: %f\n');
            offset = fscanf(file, 'Offset: %f');
            dimension = 1;
            unit = Unit(name, type, symbol, multiplier, offset, dimension);
            fclose(file);
        end
    end
end
