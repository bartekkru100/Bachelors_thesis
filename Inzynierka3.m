% Author: Bartosz Kruszona

% Variable naming guide:
% I named most of the variables after the physical properties they
% represent, most should be self-explanatory, I also added suffixes to
% signify at what station the property is being calculated or in case of
% some variables the special function it performs.

% List of suffixes:
% * at    - Atmosphere (the environment).
% * 0     - Stagnation state.
% * in    - Injection.
% * cc    - Combustion chamber after combustion.
% * th    - Throat.
% * ex    - Exit.
% * calc  - Temporary property recalculated using a different method to a
%           previously given value. Used to check for convergence.
% * old   - Value of the propery from a previous iteration. Used to check
%           for convergence.
% * min   - Lower boundary in bisection method.
% * max   - Upper boundary in bisection method.
% * 1     - Before the process. Used in a function.
% * 2     - After the process. Used in a function.

clc, clear

% Constraints
geometry_defined = false;
expansion_ratio = 80;
contraction_ratio = 5;
%{
d_cc = 0.5;
d_th  = 0.11;
d_ex = 0.11 * 2 ^ 4;
A_cc = A_circle(d_cc);
A_th = A_circle(d_th);
A_ex = A_circle(d_ex);
%}
heat_power_cc = 0;
thrust = 1.81e6;

% Atmospheric conditions
atmo = GRI30;
temperature_at = 293;
pressure_at = 1e5;

% Initial Injection
prop = Solution('CRECK.yaml');
%prop = GRI30;
temperature_in = 293;
pressure_in = 300e5;
Injection(prop, temperature_in, pressure_in)

% Initial Combustion
Combustion(prop, heat_power_cc);
pressure_cc = pressure_in;
temperature_cc = temperature(prop);
k_cc = k(prop);


% Initial stagnation state definition
pressure_0 = pressure_cc;
temperature_0 = temperature_cc;
density_0 = density(prop);
k_0 = k_cc;
entropy_0 = entropy_mass(prop);
enthalpy_0 = enthalpy_mass(prop);
a_0 = soundspeed(prop);

% Throat conditions and re-evaluation of initial and combustion states

mass_rate_area_old = 0;
for i = 1 : 20
    pressure_max = pressure_0 * (2 / (k_0 + 1)) ^ ((k_0) / (k_0 - 1)) * 1.5;
    pressure_min = pressure_max / 1.5 ^ 2;
    for j = 1 : 20
        pressure_th = (pressure_min + pressure_max) / 2;
        set(prop, 'P', pressure_th, 'S', entropy_0);
        temperature_th = temperature(prop);
        density_th = density(prop);
        enthalpy_th = enthalpy_mass(prop);
        k_th = k(prop);
        a_th = soundspeed(prop);
        velocity_th = sqrt(2 * (enthalpy_0 - enthalpy_th));
        error_V = (velocity_th - a_th) / velocity_th;
        if abs(error_V) < 1e-6 % Checking for convergence
            break;
        end
        if error_V > 0 % Range redefinition
            pressure_min = pressure_th;
        else
            pressure_max = pressure_th;
        end
    end
    % Reevaluate post-combustion state and stagnation state if mass flow
    % rate convergence (steady state condition) is not met.

    mass_rate_area = velocity_th * density_th;
    error_M = (mass_rate_area_old - mass_rate_area) / mass_rate_area_old;
    if abs(error_M) > 1e-6
        Injection(prop, temperature_in, pressure_in);
        Combustion(prop, heat_power_cc, mass_rate_area);
        pressure_cc = pressure_in;
        temperature_cc = temperature(prop);
        k_cc = k(prop);
        density_cc = density(prop);
        enthalpy_cc = enthalpy_mass(prop);
        velocity_cc = mass_rate_area / (contraction_ratio * density(prop));
        enthalpy_0 = enthalpy_mass(prop) - velocity_cc ^ 2 / 2;
        setEnthalpyIsentropic(prop, enthalpy_0);
        temperature_0 = temperature(prop);
        pressure_0 = pressure(prop);
        density_0 = density(prop);
        k_0 = k(prop);
        entropy_0 = entropy_mass(prop);
        a_0 = soundspeed(prop);
        velocity_cc = sqrt(2 * (enthalpy_0 - enthalpy_cc));
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
for i = 1 : 30
    temperature_min = 0;
    temperature_max = temperature_th;
    for j = 1 : 20
        temperature_ex = (temperature_min + temperature_max) / 2; % Bisection midpoint
        % Cantera hasn't implemented a function to change temperature at
        % constant entropy. I made my own function to faciliate that.
        setTemperatureIsentropic(prop, temperature_ex);
        pressure_ex = pressure(prop);
        density_ex = density(prop);
        k_ex1 = k(prop);
        enthalpy_ex = enthalpy_mass(prop);
        velocity_ex = sqrt(2 * (enthalpy_0 - enthalpy_ex));
        %mass_rate_calc = A_th * expansion_ratio * velocity_ex * density_ex;
        mass_rate_area_calc = expansion_ratio * velocity_ex * density_ex;
        %error_M = (mass_rate_calc - mass_rate) / mass_rate_calc;
        error_M = (mass_rate_area_calc - mass_rate_area) / mass_rate_area_calc;
        if abs(error_M) < 1e-6  % Checking for convergence
            break;
        end
        if error_M < 0 % Range redefinition
            temperature_min = temperature_ex;
        else
            temperature_max = temperature_ex;
        end
    end

    % WIP checking for expansion ratio convergence, may not be necessary,
    % as it seems that calculating it using mass flow rate, exhaust
    % velocity and density seems to give a result equal to the given
    % expansion ratio. On the other hand expansion ratio calculated using
    % the mach number relation differs slightly, but that seems to be the
    % result of k changing with temperature.
    expansion_ratio_calc = mass_rate_area / (density_ex * velocity_ex);
    mach_number_ex = velocity_ex / soundspeed(prop);
    k_ex = (k_ex1 + k_th) / 2;
    pressure_ex1 = pressure_0 * (1 + (k_ex - 1) / 2 * (mach_number_ex) ^ 2) ^ (- k_ex / (k_ex - 1));
    expansion_ratio_calc1 = ((1 + (k_ex - 1) / 2 * (mach_number_ex) ^ 2) / ((k_ex + 1) / 2)) ^ ((k_ex + 1) / (2 * (k_ex - 1))) / (mach_number_ex);
    if abs(error_M) < 1e-6
        if abs(error_M) < 1e-6
            break;
        end
    end

end
energy_ex = enthalpy_ex + velocity_ex ^ 2 / 2;
specific_impulse = (velocity_ex + (pressure_ex - 0) * expansion_ratio / mass_rate_area)% / 9.81



% Functions:

% Injection
% Cantera sets the injection composition, pressure and temperature
function Injection(prop, temperature_in, pressure_in)
%set(prop, 'Y', 'H2:1,O2:6.03', 'P', pressure_in, 'T', temperature_in);
%set(prop, 'X', 'CH4:280,O2:510', 'P', pressure_in, 'T', temperature_in);
set(prop, 'X', 'IC16H34:1', 'P', pressure_in, 'T', temperature_in);
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

% Calculates the mach number
function Mach = Mach(gas, vel)
Mach = vel/soundspeed(gas);
end

% Calculates an area of a circle
function A_circle = A_circle(d_circle)
A_circle = d_circle^2/4;
end

% Calculates a diameter of a circle
function d_circle = d_circle(A_circle)
d_circle = sqrt(A_circle) * 2;
end

% Universal gas constant
function R_universal = R_universal
R_universal = 8314.46261815324;
end

% Calculates specific gas constant
function R_specific = R_specific(gas)
R_specific = R_universal / meanMolecularWeight(gas);
end

% Calculates the specific heat ratio
function k = k(gas)
%k = 1.2;
k = cp_mass(gas)/cv_mass(gas);
end

% Sets temperature, keeping entropy constant
% Current state of the gas and target temperature of the gas are given.
% Pressure is first estimated assuming isentropic process. The initial
% estimation is inaccurate, because the specific heat ratio (k) changes with
% temperature. Bisection method is used to find an intermediate k that
% would keep the entropy constant.
function setTemperatureIsentropic(gas, temperature_2)
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
entropy_1 = entropy_mass(gas);
density_1 = density(gas);
enthalpy_1 = enthalpy_mass(gas);
pressure_1 = pressure(gas);
delta_enthalpy = enthalpy_2 - enthalpy_1;
temperature_2 = temperature(gas) + delta_enthalpy / cp_mass(gas); % dh = C_p * dT => T_2 ~= T_1 + delta_h / C_p
setTemperatureIsentropic(gas, temperature_2);
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