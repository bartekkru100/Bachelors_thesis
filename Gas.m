classdef Gas < handle

    % Property names are the same as Cantera's Solution class method names
    % expcept the _mass suffix. State class properties always refer to the
    % mass variant. Velocity property was added for convenience.

    properties (Access = public)
        solution;
    end

    properties (Access = private)
        source;
        vel;
        is_massFlowDimenstionless;
    end

    methods (Access = public)

        % Constructor takes a previously defined Solution or State Object
        % and copies its properties

        function gas = Gas(source, velocity)
            % Path for default constructor
            if nargin == 0
                gas.solution = GRI30;
                gas.source = 'GRI30.yaml';
                gas.vel = 0;

                % Path for a Solution class input
            elseif class(source) == "char"
                gas.solution = Solution(source);
                gas.source = source;
                if nargin == 2
                    gas.vel = velocity;
                else
                    gas.vel = 0;
                end

                % Path for a State class input
            elseif class(source) == "Gas"
                gas.solution = Solution(source.source);
                set(gas.solution, 'Y', source.massFractions, 'T', source.temperature, 'P', source.pressure);
                gas.source = source.source;
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

        function s_stagnation = stagnation(gas)
            import Gas.*
            state_1 = State(gas);
            s_stagnation = State(setvelocityisentropic(gas, 0));
            setstate(gas, state_1);
        end
    end

    methods (Static)

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
            Tolerance = 1e-12;
            error_S = Tolerance * 10;
            error_S_old = error_S;
            error_S_abs = 0;

            % Initial estimation of pressure after the process.

            n = 0;
            while abs(error_S) > Tolerance
                if abs(error_S_old / error_S) < 10
                    n = n + 1;
                    if n > 2
                        break;
                    end
                else
                    n = 0;
                end
                state_2_old = state_2; % Storing the state from the first estimation for later
                error_S_old = error_S;

                error_S_abs = error_S_abs + (state_1.entropy - gas.entropy);
                state_2.pressure = state_1.pressure * ((temperature_2 / state_1.temperature) ^ (state_1.cp / state_1.R_specific) / ...
                    exp(error_S_abs / state_1.R_specific));

                set(gas.solution, 'T', temperature_2, 'P', state_2.pressure);
                state_2 = State(gas);

                error_S = (state_1.entropy - gas.entropy) / state_1.entropy; % Checking for convergence
            end

            % If the estimation was not within tolerance, false position
            % method is used
            if abs(error_S) > Tolerance
                % Choosing the brackets based on the diference between
                % previous estimations
                error_P = state_2_old.pressure / state_2.pressure;

                pressure_a = state_2.pressure * error_P;
                set(gas.solution, 'T', temperature_2, 'P', pressure_a);
                error_S_a = (state_1.entropy - gas.entropy) / state_1.entropy;

                pressure_b = state_2.pressure / error_P;
                set(gas.solution, 'T', temperature_2, 'P', pressure_b);
                error_S_b = (state_1.entropy - gas.entropy) / state_1.entropy;

                m = 0;
                while 1
                    m = m + 1;
                    if m > 100
                        disp("convergence failed")
                        break;
                    end
                    state_2.pressure = (pressure_a * error_S_b - pressure_b * error_S_a) / (error_S_b - error_S_a);
                    set(gas.solution, 'T', temperature_2, 'P', state_2.pressure);
                    state_2 = State(gas);
                    error_S = (state_1.entropy - gas.entropy) / state_1.entropy; % Checking for convergence
                    if abs(error_S) < Tolerance
                        break;
                    end
                    if sign(error_S) == sign(error_S_a) % Range redefinition
                        pressure_a = state_2.pressure;
                        error_S_a = error_S;
                    else
                        pressure_b = state_2.pressure;
                        error_S_b = error_S;
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
            state_2 = State(gas);
            state_2_old = state_2;

            % Error declarations
            Tolerance = 1e-12;
            error_H = Tolerance * 10;
            error_H_old = error_H;
            error_H_abs = 0;

            % Initial estimation of pressure after the process.
            delta_enthalpy = enthalpy_2 - state_1.enthalpy;
            delta_enthalpy_calc = delta_enthalpy;

            n = 0;
            while abs(error_H) > Tolerance
                if abs(error_H_old / error_H) < 2
                    n = n + 1;
                    if n > 3
                        break;
                    end
                else
                    n = 0;
                end
                state_2_old = state_2; % Storing the state from the first estimation for later
                error_H_old = error_H;

                error_H_abs = error_H_abs + (delta_enthalpy - delta_enthalpy_calc);
                state_2.temperature = state_1.temperature + (delta_enthalpy + error_H_abs) / state_1.cp;
                state_2.pressure = state_1.pressure * (state_2.temperature / state_1.temperature) ^ (state_1.k / (state_1.k - 1));

                set(gas.solution, 'S', state_1.entropy, 'P', state_2.pressure);
                state_2 = State(gas);
                delta_enthalpy_calc = state_2.enthalpy - state_1.enthalpy;

                error_H = (enthalpy_2 - state_2.enthalpy) / enthalpy_2;
            end

            % If the estimation was not within tolerance, false position
            % method is used
            if abs(error_H) > Tolerance
                % Choosing the brackets based on the diference between
                % previous estimations
                error_P_abs = abs(state_2_old.pressure - state_2.pressure);

                pressure_a = state_2.pressure + error_P_abs;
                set(gas.solution, 'S', state_1.entropy, 'P', pressure_a);
                error_H_a = (enthalpy_2 - gas.enthalpy) / enthalpy_2;

                pressure_b = state_2.pressure - error_P_abs;
                set(gas.solution, 'S', state_1.entropy, 'P', pressure_b);
                error_H_b = (enthalpy_2 - gas.enthalpy) / enthalpy_2;

                m = 0;
                while 1
                    m = m + 1;
                    if m > 100
                        disp("convergence failed")
                        break;
                    end
                    state_2.pressure = (pressure_a * error_H_b - pressure_b * error_H_a) / (error_H_b - error_H_a);
                    set(gas.solution, 'S', state_1.entropy, 'P', state_2.pressure);
                    state_2 = State(gas);
                    error_H = (enthalpy_2 - state_2.enthalpy) / enthalpy_2; % Checking for convergence
                    if abs(error_H) < Tolerance
                        break;
                    end
                    if sign(error_H) == sign(error_H_a) % Range redefinition
                        pressure_a = state_2.pressure;
                        error_H_a = error_H;
                    else
                        pressure_b = state_2.pressure;
                        error_H_b = error_H;
                    end
                end
            end

            gas.vel = sqrt(2 * (state_1.totEnergy - gas.enthalpy));
        end

        function gas = setPressureIsentropic(gas, pressure_2)
            import Gas.*
            state1 = State(gas);
            set(gas.solution, 'P', pressure_2, 'S', gas.entropy);
            gas.vel = sqrt(2 * (state1.totEnergy - gas.enthalpy));
        end
    end
end