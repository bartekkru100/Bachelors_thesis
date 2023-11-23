classdef State

    % Property names are the same as Cantera's Solution class method names
    % expcept the _mass suffix. State class properties always refer to the
    % mass variant. Velocity property was added for convenience.

    properties (Access = public)
        temperature;
        pressure;
        density;
        meanMolecularWeight;
        enthalpy;
        intEnergy;
        entropy;
        gibbs;
        cp;
        cv;
        soundspeed;
        massFractions;
        molecularWeights;
        velocity;

    end
    methods (Access = public)

        % Constructor takes a previously defined Solution or State Object
        % and copies its properties

        function state = State(gas, velocity)
            
            % Path for default constructor
            if nargin == 0
                state.temperature = 0;
                state.pressure = 0;
                state.density = 0;
                state.meanMolecularWeight = 0;
                state.enthalpy = 0;
                state.intEnergy = 0;
                state.entropy = 0;
                state.gibbs = 0;
                state.cp = 0;
                state.cv = 0;
                state.soundspeed = 0;
                state.massFractions = 0;
                state.molecularWeights = 0;
                state.velocity = 0;

            % Path for a Solution class input
            elseif class(gas) == "Solution"
                state.temperature = temperature(gas);
                state.pressure = pressure(gas);
                state.density = density(gas);
                state.meanMolecularWeight = meanMolecularWeight(gas);
                state.enthalpy = enthalpy_mass(gas);
                state.intEnergy = intEnergy_mass(gas);
                state.entropy = entropy_mass(gas);
                state.gibbs = gibbs_mass(gas);
                state.cp = cp_mass(gas);
                state.cv = cv_mass(gas);
                state.soundspeed = soundspeed(gas);
                state.massFractions = massFractions(gas);
                state.molecularWeights = molecularWeights(gas);
                if nargin == 2
                    state.velocity = velocity;
                else
                    state.velocity = 0;
                end

            % Path for a State class input
            elseif class(gas) == "State"
                state.temperature = gas.temperature;
                state.pressure = gas.pressure;
                state.density = gas.density;
                state.meanMolecularWeight = gas.meanMolecularWeight;
                state.enthalpy = gas.enthalpy;
                state.intEnergy = gas.intEnergy;
                state.entropy = gas.entropy;
                state.gibbs = gas.gibbs;
                state.cp = gas.cp;
                state.cv = gas.cv;
                state.soundspeed = gas.soundspeed;
                state.massFractions = gas.massFractions;
                state.molecularWeights = gas.molecularWeights;
                state.velocity = gas.velocity;
            else
                throw("Invalid type.");
            end
        end
    end

    methods (Static)

        % Copies the values of either State or Solution "state" into a
        % Solution "gas"
        function setGasState(gas, state)
            import State.*
            if class(gas) == "Solution"
                if class(state) == "State"
                    set(gas, 'Y', state.massFractions, 'T', state.temperature,'P', state.pressure);
                elseif class(state) == "Solution"
                    setGasState(gas, State(state));
                else
                    throw("Invalid type of the second argument.");
                end
            else
                throw("Invalid type of the first argument.");
            end
        end

        % Calculates the specific heat ratio
        function k = k(gas)
            if class(gas) == "Solution"
                k = cp_mass(gas)/cv_mass(gas);
            elseif class(gas) == "State"
                k = gas.cp/gas.cv;
            else
                throw("Invalid type");
            end
        end

        % Calculates the mach number
        function Mach = Mach(gas, vel)
            if class(gas) == "Solution"
                Mach = vel/soundspeed(gas);
            elseif class(gas) == "State"
                Mach = vel/gas.soundspeed;
            else
                throw("Invalid type");
            end
        end

        % Universal gas constant
        function R_universal = R_universal
            R_universal = 8314.46261815324;
        end

        % Calculates specific gas constant
        function R_specific = R_specific(gas)
            import State.*
            if class(gas) == "Solution"
                R_specific = R_universal / meanMolecularWeight(gas);
            elseif class(gas) == "State"
                R_specific = R_universal / gas.meanMolecularWeight;
            else
                throw("Invalid type");
            end
        end

        % Sets temperature, keeping entropy constant
        % Current state of the gas and target temperature of the gas are given.
        % Pressure is first estimated assuming isentropic process. The initial
        % estimation is inaccurate, because the specific heat ratio (k) changes with
        % temperature. Bisection method is used to find an intermediate k that
        % would keep the entropy constant.
        function setTemperatureIsentropic(gas, temperature_2)
            import State.*
            entropy_1 = entropy_mass(gas);
            k_1 = k(gas);
            temperature_1 = temperature(gas);
            pressure_1 = pressure(gas);
            pressure_2 = pressure_1 * (temperature_2 / temperature_1) ^ (k_1 / (k_1 - 1)); % Initial pressure estimate
            set(gas, 'T', temperature_2, 'P', pressure_2);
            k_2 = k(gas);

            % While k generally drops with temperature, I used weighed average to
            % calculate k, the weights being the variable in the bisection method to
            % make convergence less likely to fail.
            weight_max = 2;
            weight_min = 0;
            i = 0;
            while 1
                i = i + 1;
                if i > 100
                    break;
                end
                weight = (weight_max + weight_min) / 2; % Bisection midpoint
                k_12 = (k_1 * (2 - weight) + k(gas) * (weight) ) / 2;
                pressure_2 = pressure_1 * (temperature_2 / temperature_1) ^ (k_12 / (k_12 - 1));
                set(gas, 'T', temperature_2, 'P', pressure_2);
                k_2 = k(gas);
                entropy_2 = entropy_mass(gas);
                error_S = (entropy_2 - entropy_1) / entropy_2; % Checking for convergence
                if abs(error_S) < 1e-6
                    break;
                end
                if error_S > 0 % Range redefinition
                    weight_min = weight;
                else
                    weight_max = weight;
                end
            end
        end

        % Sets ethalpy, keeping entropy constant
        % Works very similar to setTemperatureIsentropic, the function first
        % estimates the value of temperature after the process using the relation:
        % dh = C_p * dT, then using a relation dh = dP / rho and bisection method
        % to find an intermediate density that would keep entropy constant.
        function setEnthalpyIsentropic(gas, enthalpy_2)
            import State.*
            entropy_1 = entropy_mass(gas);
            density_1 = density(gas);
            enthalpy_1 = enthalpy_mass(gas);
            pressure_1 = pressure(gas);
            delta_enthalpy = enthalpy_2 - enthalpy_1;
            temperature_2 = temperature(gas) + delta_enthalpy / cp_mass(gas); % dh = C_p * dT => T_2 ~= T_1 + delta_h / C_p
            State.setTemperatureIsentropic(gas, temperature_2);
            density_2 = density(gas);
            weight_max = 2;
            weight_min = 0;
            i = 0;
            while 1
                i = i + 1;
                if i > 100
                    break;
                end
                weight = (weight_max + weight_min) / 2; % Bisection midpoint
                density_12 = (density_1 * (2 - weight) + density_2 * (weight) ) / 2;
                pressure_2 = pressure_1 + delta_enthalpy * density_12; % dh = dP / rho => P_2 ~= P_1 + delta_h * rho
                set(gas, 'H', enthalpy_2, 'P', pressure_2);
                density_2 = density(gas);
                entropy_2 = entropy_mass(gas);
                error_S = (entropy_2 - entropy_1) / entropy_2; % Checking for convergence
                if abs(error_S) < 1e-6
                    break;
                end
                if error_S > 0 % Range redefinition
                    weight_min = weight;
                else
                    weight_max = weight;
                end
            end
        end
    end
end