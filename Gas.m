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
            value = R_universal / meanMolecularWeight(gas.solution);;
        end

        function s_stagnation = stagnation(gas)
            import Gas.*
            s_stagnation = State(setVelocityIsentropic(gas, 0));
        end
    end

    methods (Static)

        function R_universal = R_universal()
            R_universal = 8314.46261815324;
        end

        function gas = setState(gas, property1, value1, property2, value2, property3, value3, property4, value4)
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

        function gas = setVelocityIsentropic(gas, newVelocity)
            import Gas.*
            import State.*
            newEnthalpy = gas.totEnergy - newVelocity ^ 2 / 2;
            setEnthalpyIsentropic(gas, newEnthalpy);
        end

        % Sets temperature, keeping entropy constant
        % Current state of the gas and target temperature of the gas are given.
        % Pressure is first estimated assuming isentropic process. The initial
        % estimation is inaccurate, because the specific heat ratio (k) changes with
        % temperature. Bisection method is used to find an intermediate k that
        % would keep the entropy constant.
        function gas = setTemperatureIsentropic(gas, temperature_2)
            import Gas.*
            state1 = State(gas);
            state2 = State(gas);
            state2.pressure = state1.pressure * (temperature_2 / state1.temperature) ^ (state1.k / (state1.k - 1)); % Initial pressure estimate
            set(gas.solution, 'T', temperature_2, 'P', state2.pressure);
            state2 = State(gas);

            % While k generally drops with temperature, I used weighed average to
            % calculate k, the weights being the variable in the bisection method to
            % make convergence less likely to fail.
            weight_max = 2;
            weight_min = 0;
            n = 0;
            while 1
                n = n + 1;
                if n > 100
                    break;
                end
                weight = (weight_max + weight_min) / 2; % Bisection midpoint
                k_12 = (state1.k * (2 - weight) + state2.k * (weight) ) / 2;
                state2.pressure = state1.pressure * (temperature_2 / state1.temperature) ^ (k_12 / (k_12 - 1));
                set(gas.solution, 'T', temperature_2, 'P', state2.pressure);
                state2 = State(gas);
                error_S = (state1.entropy - gas.entropy) / state1.entropy; % Checking for convergence
                if abs(error_S) < 1e-6
                    break;
                end
                if error_S < 0 % Range redefinition
                    weight_min = weight;
                else
                    weight_max = weight;
                end
            end
            gas.vel = sqrt(2 * (state1.totEnergy - gas.enthalpy));
        end

        % Sets ethalpy, keeping entropy constant
        % Works very similar to setTemperatureIsentropic, the function first
        % estimates the value of temperature after the process using the relation:
        % dh = C_p * dT, then using a relation dh = dP / rho and bisection method
        % to find an intermediate density that would keep entropy constant.
        function gas = setEnthalpyIsentropic(gas, enthalpy_2)
            import Gas.*
            import State.*
            state1 = State(gas);
            delta_enthalpy = enthalpy_2 - state1.enthalpy;
            state2.temperature = state1.temperature() + delta_enthalpy / state1.cp; % dh = C_p * dT => T_2 ~= T_1 + delta_h / C_p
            setTemperatureIsentropic(gas, state2.temperature);
            state2 = State(gas);
            weight_max = 2;
            weight_min = 0;
            n = 0;
            while 1
                n = n + 1;
                if n > 100
                    break;
                end
                weight = (weight_max + weight_min) / 2; % Bisection midpoint
                density_12 = (state1.density * (2 - weight) + state2.density * (weight) ) / 2;
                state2.pressure = state1.pressure + delta_enthalpy * density_12; % dh = dP / rho => P_2 ~= P_1 + delta_h * rho
                set(gas.solution, 'H', enthalpy_2, 'P', state2.pressure);
                state2 = State(gas);
                error_S = (state1.entropy - state2.entropy) / state1.entropy; % Checking for convergence
                if abs(error_S) < 1e-6
                    break;
                end
                if error_S < 0 % Range redefinition
                    weight_min = weight;
                else
                    weight_max = weight;
                end
            end
            gas.vel = sqrt(2 * (state1.totEnergy - gas.enthalpy));
        end

        function gas = setPressureIsentropic(gas, pressure_2)
            import Gas.*
            state1 = State(gas);
            set(gas.solution, 'P', pressure_2, 'S', gas.entropy);
            gas.vel = sqrt(2 * (state1.totEnergy - gas.enthalpy));
        end
    end
end