classdef RocketUtils
    methods (Static)
        % Injection
        % Cantera sets the injection composition, pressure and temperature
        function Injection(prop, state)
            import State.*
            setGasState(prop, state);
        end

        % Combustion
        % Cantera solves for chemical equilibrium
        % For non-chemical engines, the heat is added to the gas at constant
        % pressure
        function Combustion(prop, heat_power, mass_rate_area)
            if nargin == 3
                temperature_delta = heat_power / mass_rate_area * cp_mass(prop);
                setTemperature(prop, temperature(prop) + temperature_delta);
            end
            equilibrate(prop, 'HP');
        end
    end
end