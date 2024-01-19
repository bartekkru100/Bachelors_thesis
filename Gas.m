classdef Gas < handle

    % This class is intended to work as a shell around Cantera's own
    % Solution class, to make working with it easier for isentropic flow.
    % Additionally it can store phase diagrams. Property names are the same
    % as Cantera's Solution class method names except the _mass suffix
    % since state class properties always refer to the mass variant.

%==========================================================================

    % Object fields

    properties (Access = public)
        solution; % Object of Cantera's Solution class
        phaseDiagrams; % A phase diagram object
    end

    properties (Access = private)
        source; % Path of the source file for the solution object and phase diagram
        vel; % Velocity
    end

    methods (Access = public)

%==========================================================================

        % Constructor takes file paths for sources or a prevously defined
        % Gas object as a template.

        function gas = Gas(source)
            
            if nargin == 0 % default constructor
                source = 'GRI30';
            elseif class(source) == "Gas"
                gas.source = source.src;
            else
                gas.source = source;
            end
            gas.source = source;
            
            % The constructor looks for the phases folder inside the main
            % directory.
            sourcePath = ['Phases\', gas.source];
            gas.solution = Solution([sourcePath, '\', gas.source, '.yaml']);

            % It looks further for phase diagrams.
            sourceFolder = dir(sourcePath);
            nSubFolders = size(sourceFolder, 1);
            gas.phaseDiagrams = PhaseDiagram.empty;
            nDiagrams = 0;
            for i = 3 : nSubFolders
                pathFolder = [sourceFolder(i).folder, '\' , sourceFolder(i).name];
                if isfolder(pathFolder)
                    nDiagrams = nDiagrams + 1;
                    gas.phaseDiagrams(nDiagrams) = PhaseDiagram(pathFolder);
                end
            end

            if class(source) == "Gas"
                setstate(gas, source);
            else
                gas.vel = 0;
            end
        end

%==========================================================================

% Functions that return gas properties
% Names should be self explanatory, if not then please check the
% Cantera documentation page.

        function value = temperature(gas)
            value = temperature(gas.solution);
        end
        
        function value = pressure(gas)
            value = pressure(gas.solution);
        end

        function value = density(gas)
            value = density(gas.solution);
        end

        function value = meanMolecularWeight(gas)
            value = meanMolecularWeight(gas.solution);
        end

        function value = enthalpy(gas)
            value = enthalpy_mass(gas.solution);
        end

        function value = intEnergy(gas)
            value = intEnergy_mass(gas.solution);
        end

        function value = entropy(gas)
            value = entropy_mass(gas.solution);
        end

        function value = totEnergy(gas)
            value = enthalpy_mass(gas.solution) + gas.vel ^ 2 / 2;
        end

        function value = gibbs(gas)
            value = gibbs_mass(gas.solution);
        end

        function value = soundspeed(gas)
            value = soundspeed(gas.solution);
        end

        function value = massFractions(gas)
            value = massFractions(gas.solution);
        end

        function value = molecularWeights(gas)
            value = molecularWeights(gas.solution);
        end

        function value = velocity(gas)
            value = gas.vel;
        end

        function value = k(gas)
            value = cp_mass(gas.solution)/cv_mass(gas.solution);
        end

        function value = cp(gas)
            value = cp_mass(gas.solution);
        end

        function value = cv(gas)
            value = cv_mass(gas.solution);
        end

        function value = Mach(gas)
            value = gas.vel/soundspeed(gas.solution);
        end

        function value = R_specific(gas)
            import Gas.*
            value = R_universal / meanMolecularWeight(gas.solution);
        end

        function value = massFlowFlux(gas)
            import Gas.*
            value = gas.vel * density(gas.solution);
        end

        function value = src(gas)
            value = gas.source;
        end

%==========================================================================

    % This checks the state of matter of a given species on the phase
    % diagram

        function hasCondensation = hasCondensation(gas)
            import PhaseDiagram.*

            hasCondensation = false;
            for phaseDiagram = gas.phaseDiagrams
                if massFraction(gas.solution, phaseDiagram.speciesName) > 1e-6
                    if phaseDiagram.checkstateofmatter(gas) ~= 1
                        hasCondensation = true;
                        return;
                    end
                end
            end
        end

%==========================================================================

        % This looks for sonic conditions (v = a)

        function [s_sonic, s_stagnation] = sonic(gas, revert)
            import Gas.*

            state_1 = State(gas);
            s_stagnation = gas.stagnation;

            % Using error correction method.

            tolerance = 1e-12;
            maxIterations = 100;
            error_V_abs = -1;

            numericalMethod = ErrorCorrectedMethod("sonic_stage_1", tolerance, maxIterations);
            numericalMethod.disablewarnings;

            while 1
                error_V_abs = real(error_V_abs + (1 - gas.Mach));

                pressure = real(s_stagnation.pressure * (1  + (s_stagnation.k - 1) / 2 * (1 + error_V_abs) ^ 2) ^ ((- s_stagnation.k) / (s_stagnation.k - 1)));
                numericalMethod.findnewX(pressure);

                error_V = error_sonic(gas, pressure, state_1);

                numericalMethod.updateXY(pressure, error_V);
                if numericalMethod.checkconvergence
                    break;
                end
            end

%--------------------------------------------------------------------------

            % If the estimation was not within tolerance, Muller's method
            % is used

            if numericalMethod.hasfailed

                maxIterations = 100;

                % Point 1
                pressure(1) = numericalMethod.get.X_0;
                error_V(1) = numericalMethod.get.Y_0;

                % Point 2
                pressure(2) = (numericalMethod.get.X_0 + numericalMethod.get.X_0_old) / 2;
                error_V(2) = error_sonic(gas, pressure(2), state_1);

                % Point 3
                pressure(3) = numericalMethod.get.X_0_old;
                error_V(3) = numericalMethod.get.Y_0_old;

                numericalMethod = MullersMethod("sonic_stage_1", tolerance, maxIterations, pressure, error_V, 'max');
                numericalMethod.setX_min(0);

                while 1

                    pressure = numericalMethod.findnewX;

                    error_V = error_sonic(gas, pressure, state_1);

                    numericalMethod.updateXY(pressure, error_V);
                    if numericalMethod.checkconvergence
                        break;
                    end
                end
            end

            s_sonic = State(gas);

            if nargin == 1
            elseif revert == true
                setstate(gas, state_1);
            elseif revert == false
            else
                error("revert needs to take a logical value.");
            end
        end

%==========================================================================

        % This calculates the stagnation state

        function s_stagnation = stagnation(gas, revert)
            import Gas.*
            state_1 = State(gas);

            s_stagnation = State(setvelocityisentropic(gas, 0));

            if nargin == 1
            elseif revert == true
                setstate(gas, state_1);
            elseif revert == false
            else
                error("revert needs to take a logical value.");
            end
        end

%==========================================================================

        % This calculates the gas state for maximum possible expansion

        function s_maxVelocity = maxvelocity(gas)
            import Gas.*
            state_1 = State(gas);
            pressure = 1;
            while 1
                try
                    setpressureisentropic(gas, pressure);
                    s_maxVelocity = State(gas);
                    pressure = pressure / 10;
                catch
                    break;
                end
            end

            if nargin == 1
                setstate(gas, s_maxVelocity);
            elseif revert == true
                setstate(gas, state_1);
            elseif revert == false
                setstate(gas, s_maxVelocity);
            else
                error("revert needs to take a logical value.");
            end
        end
    end

%==========================================================================

    % WIP This imports phase diagrams from files

    methods (Access = private)
        function importphasediagram(gas)
            gas.phaseDiagrams
        end
    end

    methods (Static, Access = public)

%==========================================================================

        % Universal gas constant

        function R_universal = R_universal()
            R_universal = 8314.46261815324;
        end

%==========================================================================

        % Similar functionality to Cantera's own set() function, can also take
        % a State object for convenience.

        function gas = setstate(gas, property1, value1, property2, value2, property3, value3, property4, value4)
            import State.*
            import Gas.*

            % If one of the arguments is a State or Gas object, it copies
            % the mass fractions, temperature, pressure and velocity
            % automatically

            if ismember(class(property1), ['Gas', 'State'])
                state = property1;
                set(gas.solution, 'Y', state.massFractions, 'T', state.temperature, 'P', state.pressure);
                gas.vel = state.velocity;
                return;
            else

%--------------------------------------------------------------------------

                % Inputs are converted into two arrays that store property
                % symbol (same naming as Cantera) and the value.

                n = 0;
                if nargin > 2
                    n = n + 1;
                    propertyArr(n) = {property1};
                    valueArr(n) = {value1};

                    if nargin > 4
                        n = n + 1;
                        propertyArr(n) = {property2};
                        valueArr(n) = {value2};

                        if  nargin > 6
                            n = n + 1;
                            propertyArr(n) = {property3};
                            valueArr(n) = {value3};

                            if  nargin > 8
                                n = n + 1;
                                propertyArr(n) = {property4};
                                valueArr(n) = {value4};
                            end
                        end
                    end
                end

%--------------------------------------------------------------------------

                % Since velocity is a Gas and not a Solution property it
                % needs separate code to look for it.

                if ismember('velocity', propertyArr)
                    index = find(strcmp('velocity', propertyArr));
                    gas.vel = cell2mat(valueArr(index));
                    valueArr(index) = [];
                    propertyArr(index) = [];
                    n = n - 1;
                end

                % Values are fed into Cantera's Set() function.

                if n == 1
                    set(gas.solution, cell2mat(propertyArr(1)), cell2mat(valueArr(1)));
                elseif n == 2
                    set(gas.solution, cell2mat(propertyArr(1)), cell2mat(valueArr(1)), cell2mat(propertyArr(2)), cell2mat(valueArr(2)));
                elseif n == 3
                    set(gas.solution, cell2mat(propertyArr(1)), cell2mat(valueArr(1)), cell2mat(propertyArr(2)), cell2mat(valueArr(2)), cell2mat(propertyArr(3)), cell2mat(valueArr(3)));
                elseif n == 4
                    set(gas.solution, cell2mat(propertyArr(1)), cell2mat(valueArr(1)), cell2mat(propertyArr(2)), cell2mat(valueArr(2)), cell2mat(propertyArr(3)), cell2mat(valueArr(3)), cell2mat(propertyArr(4)), cell2mat(valueArr(4)));
                end
            end
        end

%==========================================================================

        % Sets velocity at constant entropy

        function gas = setvelocityisentropic(gas, newVelocity)
            import Gas.*
            import State.*
            newEnthalpy = gas.totEnergy - newVelocity ^ 2 / 2;
            setenthalpyisentropic(gas, newEnthalpy);
        end

%==========================================================================

        % Sets temperature at constant entropy

        function gas = settemperatureisentropic(gas, temperature_2)
            import Gas.*

            state_1 = State(gas);

            % Initial estimation of pressure after the process.

            % Using error correction method.

            tolerance = 1e-12;
            maxIterations = 10;
            error_S_abs = 0;

            numericalMethod = ErrorCorrectedMethod("settemperatureisentropic_stage1", tolerance, maxIterations);
            numericalMethod.setX_min(0);
            numericalMethod.disablewarnings;

            while 1
                error_S_abs = error_S_abs + (state_1.entropy - gas.entropy);
                pressure = state_1.pressure * ((temperature_2 / state_1.temperature) ^ (state_1.cp / state_1.R_specific) / ...
                                                      exp(error_S_abs / state_1.R_specific));
                numericalMethod.findnewX(pressure);

                error_S = error_temperatureisentropic(gas, temperature_2, pressure(2), state_1);

                numericalMethod.updateXY(state_2.pressure, error_S);
                if numericalMethod.checkconvergence()
                    break;
                end
            end

%--------------------------------------------------------------------------

            % If the estimation was not within tolerance, Muller's method
            % is used

            if numericalMethod.hasfailed

                maxIterations = 100;

                % Point 1
                pressure(1) = state_1.pressure;
                error_S(1) = error_temperatureisentropic(gas, temperature_2, pressure(1), state_1);

                % Point 2
                pressure(2) = (state_1.pressure + numericalMethod.get.X_0_old) / 2;
                error_S(2) = error_temperatureisentropic(gas, temperature_2, pressure(2), state_1);

                % Point 3
                pressure(3) = numericalMethod.get.X_0_old;
                error_S(3) = error_temperatureisentropic(gas, temperature_2, pressure(3), state_1);

                numericalMethod = MullersMethod("settemperatureisentropic_stage2", tolerance, maxIterations, pressure, error_S, 'min');
                numericalMethod.setX_min(0);

                while 1

                    pressure_0 = numericalMethod.findnewX;

                    error_S_0 = error_temperatureisentropic(gas, temperature_2, pressure_0, state_1);

                    numericalMethod.updateXY(pressure_0, error_S_0);
                    if numericalMethod.checkconvergence
                        break;
                    end
                end
            end
            % Velocity is updated using the new value of enthalpy.
            gas.vel = sqrt(2 * (state_1.totEnergy - gas.enthalpy));
        end

%==========================================================================

        % Sets ethalpy at constant entropy

        function gas = setenthalpyisentropic(gas, enthalpy_2)
            import Gas.*

            state_1 = State(gas);

            % Initial estimation of pressure after the process.

            % Using error correction method.

            tolerance = 1e-12;
            maxIterations = 10;
            error_H_abs = 0;

            delta_enthalpy = enthalpy_2 - state_1.enthalpy;
            delta_enthalpy_calc = delta_enthalpy;

            numericalMethod = ErrorCorrectedMethod("setenthalpyisentropic_stage1", tolerance, maxIterations);
            numericalMethod.setX_min(0);
            numericalMethod.disablewarnings;

            while 1

                error_H_abs = error_H_abs + (delta_enthalpy - delta_enthalpy_calc);
                temperature = state_1.temperature + (delta_enthalpy + error_H_abs) / state_1.cp;
                pressure = state_1.pressure * (temperature / state_1.temperature) ^ (state_1.k / (state_1.k - 1));

                numericalMethod.findnewX(pressure);

                error_H = error_enthalpyisentropic(gas, enthalpy_2, pressure, state_1);

                numericalMethod.updateXY(pressure, error_H);
                if numericalMethod.checkconvergence
                    break;
                end
                delta_enthalpy_calc = gas.enthalpy - state_1.enthalpy;
            end

%--------------------------------------------------------------------------

            % If the estimation was not within tolerance, Muller's method
            % is used

            if numericalMethod.hasfailed

                maxIterations = 100;

                % Point 1
                pressure(1) = state_1.pressure;
                error_H(1) = error_enthalpyisentropic(gas, enthalpy_2, pressure(1), state_1);

                % Point 2
                pressure(2) = (state_1.pressure + numericalMethod.get.X_0_old) / 2;
                error_H(2) = error_enthalpyisentropic(gas, enthalpy_2, pressure(2), state_1);

                % Point 3
                pressure(3) = numericalMethod.get.X_0_old;
                error_H(3) = error_enthalpyisentropic(gas, enthalpy_2, pressure(3), state_1);

                numericalMethod = MullersMethod("setenthalpyisentropic_stage2", tolerance, maxIterations, pressure, error_H, 'min');
                numericalMethod.setX_min(0);

                while 1
                    pressure = numericalMethod.findnewX;

                    error_H = error_enthalpyisentropic(gas, enthalpy_2, pressure, state_1);
                    
                    numericalMethod.updateXY(pressure, error_H);
                    if numericalMethod.checkconvergence
                        break;
                    end
                end
            end
            gas.vel = sqrt(2 * (state_1.totEnergy - gas.enthalpy));
        end

%==========================================================================

        % Sets pressure at constant entropy, updates velocity automatically
        % calls setpressureisentropic_still() function under the hood

        function gas = setpressureisentropic(gas, pressure_2, state_1)
            import Gas.*
            
            if nargin == 2
                state_1 = State(gas);
            end
            setpressureisentropic_still(gas, pressure_2, state_1);
            gas.vel = sqrt(2 * (state_1.totEnergy - gas.enthalpy));
        end

%--------------------------------------------------------------------------
        
        % Sets pressure at constant entropy, doesn't update velocity, but
        % checks for some potential errors while trying to calculate
        % pressure

        function gas = setpressureisentropic_still(gas, pressure_2, state_1)
            import Gas.*
            
            % Cantera has some problems calculating state variables at low
            % pressures, if the change fails the first time, it's repeated
            % in smaller increments.

            try
                set(gas.solution, 'P', pressure_2, 'S', gas.entropy);
            catch
                % Immidiately catches negative and NaN values of pressure
                if pressure_2 <= 0
                    error("Pressure must be positive.")
                elseif isnan(pressure_2)
                    error("Pressure must not be NaN.")
                end
                hasError = true;
                n = 2;
                while hasError
                    hasError = false;
                    set(gas.solution, 'T', state_1.temperature, 'P', state_1.pressure);
                    try
                        % Increments are logarithmic in scale
                        for pressure = logspace(log10(state_1.pressure), log10(pressure_2), n)
                            set(gas.solution, 'P', pressure, 'S', gas.entropy);
                        end
                    catch
                        n = n * 2;
                        if n > 64
                            error("Can't calculate isentropic pressure.")
                        end
                        hasError = true;
                    end
                end
            end
        end

    end

    methods (Static, Access = private)

%==========================================================================

        % These methods are called from other functions to make code
        % cleaner

        function error_S = error_temperatureisentropic(gas, temperature_2, pressure, state_1)
            set(gas.solution, 'T', temperature_2, 'P', pressure);
            error_S = (state_1.entropy - gas.entropy) / state_1.entropy;
        end

        function error_H = error_enthalpyisentropic(gas, enthalpy_2, pressure, state_1)
            import Gas.*

            set(gas.solution, 'T', state_1.temperature, 'P', state_1.pressure);
            setpressureisentropic_still(gas, pressure, state_1);
            error_H = (enthalpy_2 - gas.enthalpy) / enthalpy_2;
        end

        function error_V = error_sonic(gas, pressure, state_1)
            import Gas.*

            setpressureisentropic(gas, pressure, state_1);
            error_V = (gas.soundspeed - gas.velocity) / gas.soundspeed;
        end

    end
end