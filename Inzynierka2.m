clc, clear

%Geometry definitions
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

heat_power = 0;
thrust = 1.81e6;

%Initial conditions
prop = GRI30;
atmo = GRI30;
temperature_in = 293;
pressure_in = 300e5;
temperature_at = 293;
pressure_at = 1e5;

%   Injection
Injection(prop, temperature_in, pressure_in)

%   Combustion
Combustion(prop, heat_power);
pressure_cc = pressure_in;
temperature_cc = temperature(prop);
k_cc = k(prop);


%   Stagnation state
pressure_0 = pressure_cc;
temperature_0 = temperature_cc;
density_0 = density(prop);
k_0 = k_cc;
entropy_0 = entropy_mass(prop);
enthalpy_0 = enthalpy_mass(prop);
a_0 = soundspeed(prop);

%   Throat
pressure_max = 1.5 * pressure_0 * (2 / (k_cc + 1)) ^ ((k_cc) / (k_cc - 1));
pressure_min = 0.5 * pressure_max / 1.5;
%mass_rate_old = 0;
mass_rate_area_old = 0;
while 1
    pressure_th = (pressure_min + pressure_max) / 2;
    set(prop, 'P', pressure_th, 'S', entropy_0);
    temperature_th = temperature(prop);
    density_th = density(prop);
    enthalpy_th = enthalpy_mass(prop);
    k_th = k(prop);
    a_th = soundspeed(prop);
    velocity_th = sqrt(2 * (enthalpy_0 - enthalpy_th));
    %mass_rate = A_th * velocity_th * density_th;
    mass_rate_area = velocity_th * density_th;
    error_V = (velocity_th - a_th) / velocity_th;
    %error_M = (mass_rate_old - mass_rate) / mass_rate_old;
    error_M = (mass_rate_area_old - mass_rate_area) / mass_rate_area_old;
    if abs(error_V) < 1e-6
        break;
        if abs(error_M) < 1e-6
            break;
        end
    end
    if error_V > 0
        pressure_min = pressure_th;
    else
        pressure_max = pressure_th;
    end

    %       Reevaluate post-combustion state and stagnation state
    
    %{
    if abs(error_M) > 1e-6
        Injection(prop, temperature_in, pressure_in);
        Combustion(prop, heat_power, mass_rate);
        velocity_cc = mass_rate / (A_cc * density(prop));
        enthalpy_0 = enthalpy_mass(prop) - velocity_cc ^ 2 / 2;
        setState_HP(prop, [enthalpy_0, pressure_0])
        temperature_0 = temperature_cc;
        density_0 = density(prop);
        k_0 = k(prop);
        entropy_0 = entropy_mass(prop);
        a_0 = soundspeed(prop);
    end
    mass_rate_old = mass_rate;
    %}
    if abs(error_M) > 1e-6
        Injection(prop, temperature_in, pressure_in);
        Combustion(prop, heat_power, mass_rate_area);
        velocity_cc = mass_rate_area / (contraction_ratio * density(prop));
        enthalpy_0 = enthalpy_mass(prop) - velocity_cc ^ 2 / 2;
        setState_HP(prop, [enthalpy_0, pressure_0])
        temperature_0 = temperature_cc;
        density_0 = density(prop);
        k_0 = k(prop);
        entropy_0 = entropy_mass(prop);
        a_0 = soundspeed(prop);
    end
    mass_rate_area_old = mass_rate_area;
end
%expansion_ratio = A_ex / A_th;

mach_number_min = 1;
mach_number_max = 2 * sqrt(expansion_ratio);
mach_number_ex = (mach_number_min + mach_number_max) / 2;

pressure_ex = pressure_0 * (1 + (k_th - 1) / 2 * mach_number_ex ^ 2) ^ (- k_th / (k_th - 1));
set(prop, 'P', pressure_ex, 'S', entropy_0);

temperature_ex = temperature(prop);
density_ex = density(prop);
k_ex = k(prop);
enthalpy_ex = enthalpy_mass(prop);
velocity_ex = mach_number_ex * soundspeed(prop);
velocity_ex = sqrt(2 * (enthalpy_0 - enthalpy_ex));

i = 0;
while 1
    if i == 100
        break;
    end
    i = i + 1
    mach_number_ex = (mach_number_min + mach_number_max) / 2;

    pressure_ex = pressure_0 * (1 + (k_ex - 1) / 2 * (mach_number_ex) ^ 2) ^ (- k_ex / (k_ex - 1));
    set(prop, 'P', pressure_ex, 'S', entropy_0);

    temperature_ex = temperature(prop);
    density_ex = density(prop);
    k_ex = (k(prop) + k_0) / 2;
    enthalpy_ex = enthalpy_mass(prop);
    velocity_ex = mach_number_ex * soundspeed(prop);
    velocity_ex = sqrt(2 * (enthalpy_0 - enthalpy_ex));

    expansion_ratio_calc = ((1 + (k_ex - 1) / 2 * (mach_number_ex) ^ 2) / ((k_ex + 1) / 2)) ^ ((k_ex + 1) / (2 * (k_ex - 1))) / (mach_number_ex);
    %mass_rate_calc = A_th * expansion_ratio * velocity_ex * density_ex;
    mass_rate_area_calc = expansion_ratio * velocity_ex * density_ex;

    error_A = (expansion_ratio_calc - expansion_ratio) / expansion_ratio_calc;
    %error_M = (mass_rate_calc - mass_rate) / mass_rate_calc;
    error_M = (mass_rate_area_calc - mass_rate_area) / mass_rate_area_calc;
    if abs(error_A) < 1e-6
        break;
        if abs(error_M) < 1e-6
            break;
        end
    end
    if error_A < 0
        mach_number_min = mach_number_ex;
    else
        mach_number_max = mach_number_ex;
    end
end
rec = pressure_ex / pressure_th;

pressure_ex / oneatm
pressure_ex / pressure_th




%Functions:

%   Injection
function Injection(prop, temperature_in, pressure_in)
    set(prop, 'X', 'CH4:280,O2:510', 'P', pressure_in, 'T', temperature_in);
end

%   Combustion
function Combustion(prop, heat_power, mass_rate_area)
    if nargin == 3
        %temperature_delta = heat_power / (mass_rate * cp_mass(prop));
        temperature_delta = heat_power / mass_rate_area * cp_mass(prop);
        setTemperature(prop, temperature(prop) + temperature_delta);
    end
    equilibrate(prop, 'HP');
end

%   Calculates the mach number
function Mach = Mach(gas, vel)
    Mach = vel/soundspeed(gas);
end

%   Calculates an area of a circle
function A_circle = A_circle(d_circle)
    A_circle = d_circle^2/4;
end

%   Calculates a diameter of a circle
function d_circle = d_circle(A_circle)
    d_circle = sqrt(A_circle) * 2;
end

%   Universal gas constant
function R_universal = R_universal
    R_universal = 8314.46261815324;
end

%   Calculates specific gas constant
function R_specific = R_specific(gas)
    R_specific = R_universal / meanMolecularWeight(gas);
end

%   Calculates the specific heat ratio
function k = k(gas)
    %k = 1.2;
    k = cp_mass(gas)/cv_mass(gas);
end