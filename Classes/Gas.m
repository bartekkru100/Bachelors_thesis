classdef Gas < handle

    % This class is intended to work as a wrapper around Cantera's own
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

    properties (Access = public)
        source; % Path of the source file for the solution object and phase diagram
        vel; % Velocity
        s_sonic; % Stored sonic state
        s_stagnation; % Stored stagnation state
        s_maxVelocity % Stored max expansion state
    end

    methods (Access = public)

%==========================================================================

        % Constructor takes file paths for sources or a prevously defined
        % Gas object as a template.

        function gas = Gas(source)
            
            if nargin == 0 % default constructor
                source = 'GRI30';
            elseif ismember(class(source), ["Gas", "State"])
                source = source.src;
            end
            gas.source = convertStringsToChars(source);
            
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

            if ismember(class(source), ["Gas", "State"])
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

        function value = areaRatio(gas)
            import Gas.*
            if ~isempty(gas.s_sonic)
                value = gas.s_sonic.massFlowFlux / gas.massFlowFlux;
            else
                value = 0;
            end
        end

        function value = speciesNames(gas)
            value = convertCharsToStrings(speciesNames(gas.solution));
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

        % This finds equilibrium state
        function s_equilibrium = equilibrium(gas, revert)
            import Gas.*

            state_1 = State(gas);

            tolerance = 1e-9;

            equilibrate(gas.solution, 'HP', 1);

            s_equilibrium = State(gas);

            if nargin == 1
            elseif revert == true
                setstate(gas, state_1);
            elseif revert == false
            else
                error("revert needs to take a logical value.");
            end
        end

%==========================================================================

        % This looks for sonic conditions (v = a)

        function s_sonic = sonic(gas, revert)
            import Gas.*

            if ~isempty(gas.s_sonic) && (abs((gas.s_sonic.entropy - gas.entropy) / gas.s_sonic.entropy) < 1e-8)
                s_sonic = gas.s_sonic;
                if nargin == 1
                elseif revert == true
                elseif revert == false
                    setstate(gas, s_sonic);
                else
                    error("revert needs to take a logical value.");
                end
                return
            end
            
            state_1 = State(gas);
            s_stagnation = gas.stagnation;

            % Using error correction method.

            tolerance = 1e-9;
            maxIterations = 10;
            Mach_offset = 0;
            error_Mach_abs = 0;

            numericalMethod = ErrorCorrectedMethod("sonic_stage_1", tolerance, maxIterations);
            numericalMethod.disablewarnings;

            while 1
                Mach_offset = real(Mach_offset + error_Mach_abs);

                pressure = real(s_stagnation.pressure * (1  + (s_stagnation.k - 1) / 2 * (1 + Mach_offset) ^ 2) ^ ((- s_stagnation.k) / (s_stagnation.k - 1)));
                pressure = numericalMethod.findnewX(pressure);

                error_Mach = error_sonic(gas, pressure, state_1);

                numericalMethod.updateXY(pressure, error_Mach);
                if numericalMethod.checkconvergence
                    break;
                end
                error_Mach_abs = (1 - gas.Mach);
            end

%--------------------------------------------------------------------------

            % If the estimation was not within tolerance, Muller's method
            % is used

            if numericalMethod.hasfailed

                maxIterations = 100;

                % Point 1
                pressure(1) = numericalMethod.get.X_0;
                error_Mach(1) = error_sonic(gas, pressure(1), state_1);

                % Point 2
                pressure(2) = (numericalMethod.get.X_0 + numericalMethod.get.X_0_old) / 2;
                error_Mach(2) = error_sonic(gas, pressure(2), state_1);

                % Point 3
                pressure(3) = numericalMethod.get.X_0_old;
                error_Mach(3) = error_sonic(gas, pressure(3), state_1);

                numericalMethod = MullersMethod("sonic_stage_2", tolerance, maxIterations, pressure, error_Mach, 'min');
                numericalMethod.setX_min(0);

                while 1
                    pressure = numericalMethod.findnewX;

                    error_Mach = error_sonic(gas, pressure, state_1);

                    numericalMethod.updateXY(pressure, error_Mach);
                    if numericalMethod.checkconvergence
                        break;
                    end
                end
            end

            s_sonic = State(gas);
            gas.s_sonic = s_sonic;

            if nargin == 1
                setstate(gas, state_1);
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
            
            if ~isempty(gas.s_stagnation) && abs((gas.s_stagnation.entropy - gas.entropy) / gas.s_stagnation.entropy) < 1e-8
                s_stagnation = gas.s_stagnation;
                if nargin == 1
                elseif revert == true
                elseif revert == false
                    setstate(gas, s_stagnation);
                else
                    error("revert needs to take a logical value.");
                end
                return
            end
            
            state_1 = State(gas);
        
            s_stagnation = State(setvelocityisentropic(gas, 0));
            gas.s_stagnation = s_stagnation;
        
            if nargin == 1
                setstate(gas, state_1);
            elseif revert == true
                setstate(gas, state_1);
            elseif revert == false
            else
                error("revert needs to take a logical value.");
            end
        end

%==========================================================================

        % This calculates the gas state for maximum possible expansion

        function s_maxVelocity = maxvelocity(gas, revert)
            import Gas.*
            
            if ~isempty(gas.s_maxVelocity) && (abs((gas.s_maxVelocity.entropy - gas.entropy) / gas.s_maxVelocity.entropy) < 1e-8)
                s_maxVelocity = gas.s_maxVelocity;
                if nargin == 1
                elseif revert == true
                elseif revert == false
                    setstate(gas, s_maxVelocity);
                else
                    error("revert needs to take a logical value.");
                end
                return
            end

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
            
            gas.s_maxVelocity = s_maxVelocity;

            if nargin == 1
                setstate(gas, state_1);
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

%==========================================================================

    methods (Static, Access = public)

        % Universal gas constant

        function R_universal = R_universal()
            R_universal = 8314.46261815324;
        end

%==========================================================================

        % Combines 2 or more Gas objects using the same source file into
        % one

        function combined = combine(combined, gasArray, ratio, mass_mole)
            import Gas.*

            if combined == lower('new')
            combined = Gas(gasArray(1));
            end
            for gas = gasArray
                if gas.src ~= combined.src
                    error("Source files for combined gasses need to match");
                end
            end

            if nargin == 4
                if mass_mole == convertCharsToStrings(lower("mass"))
                    ratio = ratio ./ arrayfun(@(gas) gas.meanMolecularWeight, gasArray);
                elseif mass_mole == convertCharsToStrings(lower("mole"))
                else
                    error("3rd argument needs to be either be ""mass"" or ""mole""");
                end
            end
            ratio_denominator = sum(ratio);

            pressure_partial = arrayfun(@(gas) gas.pressure, gasArray) .* ratio  / ratio_denominator;
            pressure_total = sum(pressure_partial);

            intEnergy_partial = arrayfun(@(gas) gas.intEnergy, gasArray) .* pressure_partial / pressure_total;
            intEnergy_total = sum(intEnergy_partial);

            massFractions_partial = arrayfun(@(gas) {gas.massFractions * gas.meanMolecularWeight}, gasArray);
            massFractions_partial = cell2mat(massFractions_partial);
            massFractions_partial = massFractions_partial .* (pressure_partial / pressure_total);
            massFractions_total = sum(massFractions_partial, 2);

            density_partial = arrayfun(@(gas) gas.density, gasArray) .* ratio  / ratio_denominator;
            density_total = sum(density_partial);
            velocity_partial = arrayfun(@(gas) (gas.density * gas.velocity), gasArray) .* ratio  / ratio_denominator;
            velocity_total = sum(velocity_partial);
            setstate(combined, 'Rho', density_total, 'P', pressure_total, 'Y', massFractions_total, 'velocity', velocity_total);
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

            % Using error correction method.

            tolerance = 1e-9;
            maxIterations = 10;
            S_offset = 0;
            error_S_abs = 0;

            numericalMethod = ErrorCorrectedMethod("settemperatureisentropic_stage1", tolerance, maxIterations);
            numericalMethod.setX_min(0);
            numericalMethod.disablewarnings;

            while 1
                S_offset = S_offset + error_S_abs;
                pressure = state_1.pressure * ((temperature_2 / state_1.temperature) ^ (state_1.cp / state_1.R_specific) / ...
                                                      exp(S_offset / state_1.R_specific));
                pressure = numericalMethod.findnewX(pressure);

                error_S = error_temperatureisentropic(gas, temperature_2, pressure, state_1);

                numericalMethod.updateXY(pressure, error_S);
                if numericalMethod.checkconvergence()
                    break;
                end
                error_S_abs = state_1.entropy - gas.entropy;
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
            delta_enthalpy = enthalpy_2 - state_1.enthalpy;

            % Using error correction method.

            tolerance = 1e-9;
            maxIterations = 10;
            error_H_abs = 0;

            numericalMethod = ErrorCorrectedMethod("setenthalpyisentropic_stage1", tolerance, maxIterations);
            numericalMethod.setX_min(0);
            numericalMethod.disablewarnings;

            while 1
                delta_enthalpy = delta_enthalpy + error_H_abs;
                temperature = state_1.temperature + (delta_enthalpy) / state_1.cp;
                pressure = state_1.pressure * (temperature / state_1.temperature) ^ (state_1.k / (state_1.k - 1));

                pressure = numericalMethod.findnewX(pressure);

                error_H = error_enthalpyisentropic(gas, enthalpy_2, pressure, state_1);

                numericalMethod.updateXY(pressure, error_H);
                if numericalMethod.checkconvergence
                    break;
                end
                error_H_abs = enthalpy_2 - gas.enthalpy;
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
        
        % This function calculates the exit state assuming a completely supersonic
        % flow in the diverging part of the nozzle. Iterating over pressure to
        % satisfy mass conservation.
    
        function gas = setarearatioisentropic(gas, areaRatio_2, sub_super)
            import Gas.*
        
            s_sonic = gas.sonic;
            s_stagnation = gas.stagnation;

            % Using Muller's method
            tolerance = 1e-9;
            iterationLimit = 100;
            if sub_super == "sub"
            % Point 1
            pressure(1) = s_stagnation.pressure * (1 - 1 / areaRatio_2 ^ 2.3) + s_sonic.pressure * (1 / areaRatio_2 ^ 2.3);
            error_M(1) = error_areaisentropic(gas, areaRatio_2, pressure(1));
        
            % Point 2
            pressure(2) = s_stagnation.pressure * (1 - 1 / areaRatio_2 ^ 2.5) + s_sonic.pressure * (1 / areaRatio_2 ^ 2.5);
            error_M(2) = error_areaisentropic(gas, areaRatio_2, pressure(2));
        
            % Point 3
            pressure(3) = s_stagnation.pressure * (1 - 1 / areaRatio_2 ^ 2.7) + s_sonic.pressure * (1 / areaRatio_2 ^ 2.7);
            error_M(3) = error_areaisentropic(gas, areaRatio_2, pressure(3));
        
            numericalMethod = MullersMethod("setsupersonicexitconditions", tolerance, iterationLimit, pressure, error_M, 'max');
            numericalMethod.setX_min_max(s_sonic.pressure, s_stagnation.pressure);
        
            else
            % Point 1
            pressure(1) = s_sonic.pressure * (1 / areaRatio_2 ^ 1.8);
            error_M(1) = error_areaisentropic(gas, areaRatio_2, pressure(1));
        
            % Point 2
            pressure(2) = s_sonic.pressure * (1 / areaRatio_2 ^ 2);
            error_M(2) = error_areaisentropic(gas, areaRatio_2, pressure(2));
        
            % Point 3
            pressure(3) = s_sonic.pressure * (1 / areaRatio_2 ^ 2.2);
            error_M(3) = error_areaisentropic(gas, areaRatio_2, pressure(3));
        
            numericalMethod = MullersMethod("setsupersonicexitconditions", tolerance, iterationLimit, pressure, error_M, 'min');
            numericalMethod.setX_min_max(0, s_sonic.pressure);
            end
            while 1
                pressure = numericalMethod.findnewX;
        
                error_M = error_areaisentropic(gas, areaRatio_2, pressure);
        
                numericalMethod.updateXY(pressure, error_M);
                if numericalMethod.checkconvergence
                    break;
                end
            end
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
            setpressureisentropic_still(gas, pressure, state_1.pressure);
            error_H = (enthalpy_2 - gas.enthalpy) / enthalpy_2;
        end

        function error_Mach = error_sonic(gas, pressure, state_1)
            import Gas.*

            setpressureisentropic(gas, pressure, state_1);
            error_Mach = 1 - gas.Mach;
        end

        function error_M = error_areaisentropic(gas, areaRatio_2, pressure)
            import Gas.*
            setstate(gas, gas.s_sonic);
            setpressureisentropic(gas, pressure, gas.s_sonic);

            massFlowFlux_calc = areaRatio_2 * gas.massFlowFlux;
            error_M = (massFlowFlux_calc - gas.s_sonic.massFlowFlux) / massFlowFlux_calc;
        end

    end
end