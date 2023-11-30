% Author: Bartosz Kruszona
import State.*

clc, clear

% Constraints
geometry_defined = false;
expansion_ratio = 80;
contraction_ratio = 5;
heat_power_area = 0;
thrust = 1;

% Atmospheric conditions
atmo = GRI30;
temperature_at = 293;
pressure_at = 0.1;


% Initial Injection
prop = Solution('R1highT.yaml');
%prop = GRI30;  
%set(prop, 'Y', 'H2:1,O2:6.03');
set(prop, 'Y', 'CH4:1,O2:3.6');
%set(prop, 'X', 'N2:1');
%set(prop, 'Y', 'POSF7688:1,O2:2.72');
s_injection = State(prop);
s_injection.pressure = 300e5;
s_injection.temperature = 293;
injection(prop, s_injection);
% Initial Combustion
combustion(prop);
s_chamber = State(prop);



% Initial stagnation state definition
s_stagnation = State(s_chamber);

% Throat conditions and re-evaluation of initial and combustion states
s_throat = State();
mass_rate_area_old = 0;
for i = 1 : 20
    pressure_max = s_stagnation.pressure * (2 / (k(s_stagnation) + 1)) ^ ((k(s_stagnation)) / (k(s_stagnation) - 1)) * 1.5;
    pressure_min = s_stagnation.pressure / 1.5 ^ 2;
    for j = 1 : 20
        s_throat.pressure = (pressure_min + pressure_max) / 2;
        set(prop, 'P', s_throat.pressure, 'S', s_stagnation.entropy);
        s_throat = State(prop);
        s_throat.velocity = sqrt(2 * (s_stagnation.enthalpy - s_throat.enthalpy));
        error_V = (s_throat.velocity - s_throat.soundspeed) / s_throat.velocity;
        if abs(error_V) < 1e-6 % Checking for convergence
            break;
        end
        if error_V > 0 % Range redefinition
            pressure_min = s_throat.pressure;
        else
            pressure_max = s_throat.pressure;
        end
    end
    % Reevaluate post-combustion state and stagnation state if mass flow
    % rate convergence (steady state condition) is not met.

    mass_rate_area = s_throat.velocity * s_throat.density;
    error_M = (mass_rate_area_old - mass_rate_area) / mass_rate_area_old;
    if abs(error_M) > 1e-6
        injection(prop, s_injection);
        combustion(prop, heat_power_area, mass_rate_area);
        s_chamber = State(prop);
        s_chamber.velocity = mass_rate_area / (contraction_ratio * s_chamber.density);
        s_stagnation.enthalpy = s_chamber.enthalpy - s_chamber.velocity ^ 2 / 2;
        setEnthalpyIsentropic(prop, s_stagnation.enthalpy);
        s_stagnation = State(prop);
        s_chamber.velocity = sqrt(2 * (s_chamber.enthalpy - s_stagnation.enthalpy));
    else
        if abs(error_V) < 1e-6  % Checking for convergence
            break;
        end
    end
    mass_rate_area_old = mass_rate_area;
end

% Bisection method, used to find the exit state. Temperature is used as a
% variable and is known to fall between absolute zero and throat
% temperature. Mass flow rate is used to checked for convergence to satisfy
% conservation of mass.
s_exit = State();
for i = 1 : 30
    temperature_min = 0;
    temperature_max = s_throat.temperature;
    for j = 1 : 20
        s_exit.temperature = (temperature_min + temperature_max) / 2; % Bisection midpoint
        % Cantera hasn't implemented a function to change temperature at
        % constant entropy. I made my own function to faciliate that.
        setTemperatureIsentropic(prop, s_exit.temperature);
        s_exit = State(prop);
        s_exit.velocity = sqrt(2 * (s_stagnation.enthalpy - s_exit.enthalpy));
        mass_rate_area_calc = expansion_ratio * s_exit.velocity * s_exit.density;
        error_M = (mass_rate_area_calc - mass_rate_area) / mass_rate_area_calc;
        if abs(error_M) < 1e-6  % Checking for convergence
            break;
        end
        if error_M < 0 % Range redefinition
            temperature_min = s_exit.temperature;
        else
            temperature_max = s_exit.temperature;
        end
    end
end
specific_impulse = (s_exit.velocity + (s_exit.pressure - pressure_at) * expansion_ratio / mass_rate_area) / 9.81
mass_rate = thrust / specific_impulse;
A_th = mass_rate / mass_rate_area;
D_th = d_circle(A_th)
heat_power = heat_power_area * A_th;

% Functions:


% Calculates an area of a circle
function A_circle = A_circle(d_circle)
A_circle = d_circle^2/4;
end

% Calculates a diameter of a circle
function d_circle = d_circle(A_circle)
d_circle = sqrt(A_circle) * 2;
end