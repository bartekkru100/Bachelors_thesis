% Author: Bartosz Kruszona
import State.*
import Gas.*

clc, clear

% Constraints
geometry_defined = false;
expansion_ratio = 36.87;
contraction_ratio = 11.7;
heat_power_area = 0;
thrust = 3.83e6;

% Atmospheric conditions
atmo = Gas(Air);
temperature_at = 293;
pressure_at = oneatm;


% Initial Injection
pressure_in = 26.7e6;
temperature_in = 100;
prop = Gas('R1highT.yaml');
%prop = Gas('GRI30.yaml');  
%setState(prop, 'Y', 'H2:1,O2:6.03');
%setState(prop, 'Y', 'CH4:1,O2:3.6');
%setState(prop, 'X', 'N2:1');
setState(prop, 'Y', 'POSF7688:1,O2: 2.36');
setState(prop, 'P', pressure_in, 'T', temperature_in);
s_injection = State(prop);
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
    pressure_max = s_stagnation.pressure * (2 / (s_stagnation.k + 1)) ^ ((s_stagnation.k) / (s_stagnation.k - 1)) * 1.5;
    pressure_min = pressure_max / 1.5 ^ 2;
    for j = 1 : 20
        s_throat.pressure = (pressure_min + pressure_max) / 2;
        setPressureIsentropic(prop, s_throat.pressure);
        s_throat = State(prop);
        error_V = (prop.velocity - prop.soundspeed) / prop.velocity;
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

    massRate_area = s_throat.velocity * s_throat.density;
    error_M = (mass_rate_area_old - massRate_area) / mass_rate_area_old;
    if abs(error_M) > 1e-6
        injection(prop, s_injection);
        combustion(prop, heat_power_area, massRate_area);
        s_chamber = State(prop);
        s_chamber.velocity = massRate_area / (contraction_ratio * s_chamber.density);
        s_stagnation.enthalpy = s_chamber.enthalpy - s_chamber.velocity ^ 2 / 2;
        setEnthalpyIsentropic(prop, s_stagnation.enthalpy);
        s_stagnation = State(prop);
        s_chamber.velocity = sqrt(2 * (s_chamber.enthalpy - s_stagnation.enthalpy));
    else
        if abs(error_V) < 1e-6  % Checking for convergence
            break;
        end
    end
    mass_rate_area_old = massRate_area;
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
        massRate_area_calc = expansion_ratio * s_exit.velocity * s_exit.density;
        error_M = (massRate_area_calc - massRate_area) / massRate_area_calc;
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
s_exit.Mach
specific_impulse = (s_exit.velocity + (s_exit.pressure - pressure_at) * expansion_ratio / massRate_area)% / 9.81
mass_rate = thrust / specific_impulse;
A_th = mass_rate / massRate_area;
A_ex = A_th * expansion_ratio;
D_th = d_circle(A_th)
D_ex = d_circle(A_ex)
heat_power = heat_power_area * A_th;

% Functions:


% Calculates an area of a circle
function A_circle = A_circle(d_circle)
A_circle = d_circle^2/4;
end

% Calculates a diameter of a circle
function d_circle = d_circle(A_circle)
d_circle = sqrt(A_circle / pi) * 2;
end