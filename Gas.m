classdef Gas < handle

    % Property names are the same as Cantera's Solution class method names
    % expcept the _mass suffix. State class properties always refer to the
    % mass variant. Velocity property was added for convenience.

    properties (Access = public)
        solution;
        phaseDiagrams;
    end

    properties (Access = private)
        solutionSource;
        phaseDiagramSource;
        vel;
        isMassFlowDimenstionless;
    end

    methods (Access = public)

        % Constructor takes a previously defined Solution or State Object
        % and copies its properties

        function gas = Gas(source, velocity)
            % Path for default constructor

            gas.phaseDiagrams = PhaseDiagram('D:\Uni\Inzynierka\Program\PhaseCurves\H2O');
            if nargin == 0
                gas.solution = GRI30;
                gas.solutionSource = 'GRI30.yaml';
                gas.vel = 0;

                % Path for a Solution class input
            elseif class(source) == "char"
                gas.solution = Solution(source);
                gas.solutionSource = source;
                if nargin == 2
                    gas.vel = velocity;
                else
                    gas.vel = 0;
                end

                % Path for a State class input
            elseif class(source) == "Gas"
                gas.solution = Solution(source.solutionSource);
                set(gas.solution, 'Y', source.massFractions, 'T', source.temperature, 'P', source.pressure);
                gas.solutionSource = source.solutionSource;
                gas.vel = source.velocity;
            end
        end

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

        function hasCondensation = hasCondensation(gas)
            import PhaseDiagram.*

            hasCondensation = false;
            for i = 1 : size(gas.phaseDiagrams, 2)
                if gas.massFractions ~= 0
                    if gas.phaseDiagrams.checkstateofmatter(gas) ~= 1
                        hasCondensation = true;
                        return;
                    end
                end
            end
        end

        function s_stagnation = stagnation(gas)
            import Gas.*
            state_1 = State(gas);
            s_stagnation = State(setvelocityisentropic(gas, 0));
            setstate(gas, state_1);
        end

        function s_macVelocity = maxvelocity(gas)
            import Gas.*
            state_1 = State(gas);
            pressure = 1;
            while 1
                try
                    setpressureisentropic(gas, pressure);
                    s_macVelocity = State(gas);
                    pressure = pressure / 10;
                catch
                    setstate(gas, state_1);
                    return;
                end
            end
        end
    end

    methods (Access = private)
        function importphasediagram(gas)
            gas.phaseDiagrams
        end
    end

    methods (Static, Access = public)

        function R_universal = R_universal()
            R_universal = 8314.46261815324;
        end

        function gas = setstate(gas, property1, value1, property2, value2, property3, value3, property4, value4)
            import State.*
            import Gas.*
            if ismember(class(property1), ['Gas', 'State'])
                state = property1;
                set(gas.solution, 'Y', state.massFractions, 'T', state.temperature, 'P', state.pressure);
                gas.vel = state.velocity;
                return;
            else
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

                if ismember('velocity', propertyArr)
                    index = find(strcmp('velocity', propertyArr));
                    gas.vel = cell2mat(valueArr(index));
                    valueArr(index) = [];
                    propertyArr(index) = [];
                    n = n - 1;
                end

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

        function gas = setvelocityisentropic(gas, newVelocity)
            import Gas.*
            import State.*
            newEnthalpy = gas.totEnergy - newVelocity ^ 2 / 2;
            setenthalpyisentropic(gas, newEnthalpy);
        end

        % Sets temperature, keeping entropy constant
        % Current state of the gas and target temperature of the gas are given.
        % Pressure is first estimated assuming isentropic process. The initial
        % estimation is inaccurate, because the specific heat ratio (k) changes with
        % temperature. Bisection method is used to find an intermediate k that
        % would keep the entropy constant.
        function gas = settemperatureisentropic(gas, temperature_2)
            import Gas.*

            % Declaration of states.
            state_1 = State(gas);
            state_2 = State(gas);
            state_2_old = state_2;

            % Error declarations
            tolerance = 1e-12;
            error_S_abs = 0;

            % Initial estimation of pressure after the process.

            numericalMethod = ErrorCorrectedMethod("settemperatureisentropic_stage1", tolerance, 100);

            while 1
                error_S_abs = error_S_abs + (state_1.entropy - gas.entropy);
                pressure = state_1.pressure * ((temperature_2 / state_1.temperature) ^ (state_1.cp / state_1.R_specific) / ...
                                                      exp(error_S_abs / state_1.R_specific));
                numericalMethod.findnewX(pressure);

                set(gas.solution, 'T', temperature_2, 'P', pressure);
                state_2 = State(gas);

                error_S = (state_1.entropy - gas.entropy) / state_1.entropy; % Checking for convergence

                numericalMethod.updateXY(state_2.pressure, error_S);
                if numericalMethod.checkconvergence()
                    break;
                end
            end

            % If the estimation was not within tolerance, Muller's position
            % method is used

            if ~numericalMethod.hasconverged
                
                % Point 1
                pressure(1) = state_2.pressure;
                error_S(1) = error_temperatureisentropic(gas, temperature_2, pressure(1), state_1);
                % Point 2
                pressure(2) = (state_2.pressure + numericalMethod.get.X_0_old) / 2;
                error_S(2) = error_temperatureisentropic(gas, temperature_2, pressure(2), state_1);

                % Point 3
                pressure(3) = numericalMethod.get.X_0_old;
                error_S(3) = error_temperatureisentropic(gas, temperature_2, pressure(3), state_1);

                numericalMethod = MullersMethod("settemperatureisentropic_stage2", tolerance, 100, pressure, error_S, 'min');

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

        % Sets ethalpy, keeping entropy constant
        % Works very similar to setTemperatureIsentropic, the function first
        % estimates the value of temperature after the process using the relation:
        % dh = C_p * dT, then using a relation dh = dP / rho and bisection method
        % to find an intermediate density that would keep entropy constant.

        function gas = setenthalpyisentropic(gas, enthalpy_2)
            import Gas.*

            % Declaration of states.
            state_1 = State(gas);
            pressure_min = 0.1;%gas.maxvelocity.pressure;

            % Error declarations
            tolerance = 1e-12;
            maxIterations = 10;
            error_H_abs = 0;

            % Initial estimation of pressure after the process.
            delta_enthalpy = enthalpy_2 - state_1.enthalpy;
            delta_enthalpy_calc = delta_enthalpy;

            numericalMethod = ErrorCorrectedMethod("setenthalpyisentropic_stage1", tolerance, 100);

            while 1

                error_H_abs = error_H_abs + (delta_enthalpy - delta_enthalpy_calc);
                temperature = state_1.temperature + (delta_enthalpy + error_H_abs) / state_1.cp;
                pressure = state_1.pressure * (temperature / state_1.temperature) ^ (state_1.k / (state_1.k - 1));

                numericalMethod.findnewX(pressure);

                setpressureisentropic_still(gas, pressure, state_1);
                delta_enthalpy_calc = gas.enthalpy - state_1.enthalpy;

                error_H = (enthalpy_2 - gas.enthalpy) / enthalpy_2;

                numericalMethod.updateXY(pressure, error_H);
                if numericalMethod.checkconvergence
                    break;
                end
            end

            % If the estimation was not within tolerance, Muller's position
            % method is used
            if ~numericalMethod.hasconverged

                pressure_min = gas.maxvelocity.pressure;

                % Point 1
                pressure(1) = pressure;
                error_H(1) = error_enthalpyisentropic(gas, enthalpy_2, pressure(1), state_1);

                % Point 2
                pressure(2) = (pressure(1) + numericalMethod.get.Y_0_old) / 2;
                error_H(2) = error_enthalpyisentropic(gas, enthalpy_2, pressure(2), state_1);

                % Point 3
                pressure(3) = numericalMethod.get.Y_0_old;
                error_H(3) = error_enthalpyisentropic(gas, enthalpy_2, pressure(3), state_1);

                numericalMethod = MullersMethod("setenthalpyisentropic_stage2", tolerance, 100, pressure, error_H, 'min');

                while 1
                    pressure = numericalMethod.findnewX;
                    if pressure <= 0
                        pressure = pressure_min;
                    end
                    error_H = error_enthalpyisentropic(gas, enthalpy_2, pressure, state_1);
                    numericalMethod.updateXY(pressure, error_H);
                    if numericalMethod.checkconvergence
                        break;
                    end
                end
            end
            gas.vel = sqrt(2 * (state_1.totEnergy - gas.enthalpy));
        end

        function gas = setpressureisentropic(gas, pressure_2, state_1)
            import Gas.*
            
            if nargin == 2
                state_1 = State(gas);
            end
            setpressureisentropic_still(gas, pressure_2, state_1);
            gas.vel = sqrt(2 * (state_1.totEnergy - gas.enthalpy));
        end

        function gas = setpressureisentropic_still(gas, pressure_2, state_1)
            import Gas.*

            try
                set(gas.solution, 'P', pressure_2, 'S', gas.entropy);
            catch
                if pressure_2 <= 0
                    error("Pressure must be positive.")
                elseif isnan(pressure_2)
                    error("Pressure must not be NaN.")
                end
                hasError = true;
                n = 2;
                while hasError
                    hasError = false;
                    setstate(gas, state_1);
                    try
                        for pressure = logspace(log10(state_1.pressure), log10(pressure_2), n)
                            set(gas.solution, 'P', pressure, 'S', gas.entropy);
                        end
                    catch
                        n = n * 2;
                        if n > 1024
                            error("Can't calculate isentropic pressure.")
                        end
                        hasError = true;
                    end
                end
            end
        end

    end

    methods (Static, Access = private)

        function error_S = error_temperatureisentropic(gas, temperature_2, pressure, state_1)
            import Gas.*

            set(gas.solution, 'T', temperature_2, 'P', pressure);
            error_S = (state_1.entropy - gas.entropy) / state_1.entropy; % Checking for convergence
        end

        function error_H = error_enthalpyisentropic(gas, enthalpy_2, pressure, state_1)
            import Gas.*

            setstate(gas, state_1);
            setpressureisentropic_still(gas, pressure, state_1);
            error_H = (enthalpy_2 - gas.enthalpy) / enthalpy_2;
        end

    end
end