% Author: Bartosz Kruszona
import State.*
import Gas.*

clc, clear

% Constraints
geometry_defined = false;
expansionRatio = 34.34;
contractionRatio = 10;
heat_power_area = 0;
thrust = 2.26e6;
Tolerance = 1e-12;

% Atmospheric conditions
atmo = Gas('air.yaml');
temperature_at = 293;
pressure_at = oneatm;
setstate(atmo, 'P', pressure_at, 'T', temperature_at);
s_atmo = State(atmo);


% Initial Injection
pressure_in = 300e5;
temperature_in = 150;
gas = Gas('R1highT.yaml');
%prop = Gas('GRI30.yaml');  
%setState(prop, 'Y', 'H2:1,O2:6.03');
setstate(gas, 'Y', 'CH4:22,O2:78');
%setstate(prop, 'X', 'N2:1');
%setstate(prop, 'Y', 'POSF7688:1,O2: 2.36');
setstate(gas, 'P', pressure_in, 'T', temperature_in);
s_injection = State(gas);
% Initial Combustion
setchamberconditions(gas, s_injection);
s_chamber = State(gas);




% Initial stagnation state definition

% Throat conditions and re-evaluation of initial and combustion states
mass_rate_area_old = 0;
for i = 1 : 20
    setthroatconditions(gas);
    s_throat = State(gas);

    % Reevaluate post-combustion state and stagnation state if mass flow
    % rate convergence (steady state condition) is not met.

    massFlowFlux = s_throat.massFlowFlux;
    error_M = (mass_rate_area_old - massFlowFlux) / mass_rate_area_old;
    if abs(error_M) < Tolerance
        break;
    end

    setinjectionconditions(gas, s_injection, contractionRatio);
    s_injection = State(gas);

    setchamberconditions(gas, s_injection, heat_power_area);
    s_chamber = State(gas);

    mass_rate_area_old = massFlowFlux;
end

s_stagnation = gas.stagnation;

% Bisection method, used to find the exit state. Temperature is used as a
% variable and is known to fall between absolute zero and throat
% temperature. Mass flow rate is used to checked for convergence to satisfy
% conservation of mass.
s_exit = State();
for i = 1 : 1
    temperature_min = 0;
    temperature_max = s_throat.temperature;
    for j = 1 : 200
        setstate(gas, s_throat);
        s_exit.temperature = (temperature_min + temperature_max) / 2; % Bisection midpoint
        % Cantera hasn't implemented a function to change temperature at
        % constant entropy. I made my own function to faciliate that.
        settemperatureisentropic(gas, s_exit.temperature);
        s_exit = State(gas);
        massRate_area_calc = expansionRatio * s_exit.velocity * s_exit.density;
        error_M = (massRate_area_calc - s_throat.massFlowFlux) / massRate_area_calc;
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
if s_stagnation.pressure < s_atmo.pressure;
    shock = shockposition(gas, s_atmo, s_exit, s_throat);
    shock * gas.density * gas.velocity;
end
s_exit = State(gas);
s_exit.Mach;
specific_impulse = (s_exit.velocity + (s_exit.pressure - pressure_at) * expansionRatio / s_throat.massFlowFlux);% / 9.81
mass_rate = thrust / specific_impulse;
A_th = mass_rate / s_throat.massFlowFlux;
A_ex = A_th * expansionRatio;
D_th = d_circle(A_th);
D_ex = d_circle(A_ex);
heat_power = heat_power_area * A_th;
disp("Velocity: " + s_injection.velocity + " Entropy:" + s_injection.entropy + " Total Energy: " + s_injection.totEnergy + " Mass Flow Rate: " + s_injection.massFlowFlux * contractionRatio)
disp("Velocity: " + s_chamber.velocity + " Entropy: " + s_chamber.entropy + " Total Energy: " + s_chamber.totEnergy + " Mass Flow Rate: " + s_chamber.massFlowFlux * contractionRatio)
disp("Velocity: " + s_throat.velocity + " Entropy: " + s_throat.entropy + " Total Energy: " + s_throat.totEnergy + " Mass Flow Rate: " + s_throat.massFlowFlux)
disp("Velocity: " + s_exit.velocity + " Entropy: " + s_exit.entropy + " Total Energy: " + s_exit.totEnergy + " Mass Flow Rate: " + s_exit.massFlowFlux * expansionRatio)

% Functions:

% Calculates an area of a circle
function A_circle = A_circle(d_circle)
A_circle = d_circle^2/4;
end

% Calculates a diameter of a circle
function d_circle = d_circle(A_circle)
d_circle = sqrt(A_circle / pi) * 2;
end