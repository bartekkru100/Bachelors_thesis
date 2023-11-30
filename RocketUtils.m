classdef RocketUtils
    methods (Static)
        % Injection
        % Cantera sets the injection composition, pressure and temperature
        function injection(prop, state)
            import State.*
            setgasstate(prop, state);
        end

        % Combustion
        % Cantera solves for chemical equilibrium
        % For non-chemical engines, the heat is added to the gas at constant
        % pressure
        function combustion(prop, heat_power, mass_rate_area)
            import State.*
            equilibrate(prop, 'HP');
            if nargin == 1
                return;
            end
            heat_delta = heat_power / mass_rate_area;
            temperature_delta = heat_delta / cp_mass(prop);
            s_preHeat = State(prop);
            setTemperature(prop, s_preHeat.temperature + temperature_delta);
            s_postHeat = State(prop);
            temperature_min = s_preHeat.temperature;
            temperature_max = s_postHeat.temperature;
            for i = 1:20
                temperature_12 = (temperature_min + temperature_max) / 2;
                set(prop, 'T', temperature_12, 'P', s_preHeat.pressure);
                cp_12 = cp_mass(prop);
                setgasstate(prop, s_preHeat);
                temperature_delta = heat_delta / cp_12;
                set(prop, 'T', s_preHeat.temperature + temperature_delta, 'P', s_preHeat.pressure);
                s_postHeat = State(prop);
                error_S = (s_postHeat.entropy - (s_preHeat.entropy + heat_delta / temperature_12)) / s_postHeat.entropy; % Checking for convergence
                if abs(error_S) < 1e-6
                    break;
                end
                if error_S < 0 % Range redefinition
                    temperature_min = temperature_12;
                else
                    temperature_max = temperature_12;
                end
            end
        end
    end
end