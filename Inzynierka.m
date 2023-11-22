clc, clear, cleanup


%Geometry definitions
d_ch = 0.5;
d_th  = 0.11;
d_ex = 0.11 * 2 ^ 1;
A_ch = A_circle(d_ch);
A_th = A_circle(d_th);
A_ex = A_circle(d_ex);

%Initial conditions
prop = GRI30;
atmo = GRI30;
temperature_in = 293;
pressure_in = 300e5;
temperature_at = 293;
pressure_at = 1e5;
temperature_cc = 3950;
pressure_cc = pressure_in;

%Initialization

%   Injection
set(prop, 'X', 'CH4:280,O2:510', 'P', pressure_in, 'T', temperature_in);
kRecord(1) = k(prop);
m_rate = 1;
velocity_in = m_rate / (A_ch * density(prop));

%   Combustion/heating
set(prop, 'P', pressure_cc, 'T', temperature_cc);
kRecord(2) = k(prop);
equilibrate(prop, 'TP');
velocity_ch = m_rate / (A_ch * density(prop));

%   Stagnation state
pressure_0 = pressure_cc;
temperature_0 = temperature_cc;
density_0 = density(prop);
k_0 = k(prop);
entropy_0 = entropy_mass(prop);
enthalpy_0 = enthalpy_mass(prop) + velocity_ch ^ 2;
prop_0 = prop;
a_0 = soundspeed(prop_0);

%   Throat
m_rate = A_th * sqrt(pressure(prop) * density(prop) * k(prop) * (1 + (k(prop) - 1) / 2) ^ ( - (k(prop) + 1) / (k(prop) - 1)));
pressure_th = pressure(prop) * (2 / (k(prop) + 1)) ^ ((k(prop)) / (k(prop) - 1));
set(prop, 'P', pressure_th, 'S', entropy_0);

k_max = k_0;
pressure_max = 2 * pressure(prop);
k_min = 2 * k(prop) - k_max;
pressure_min = pressure_max / 2;

%       bisection method for throat pressure
while 1
    pressure_th = (pressure_min + pressure_max) / 2;
    set(prop, 'P', pressure_th, 'S', entropy_0);
    m_rate_calc = A_th * soundspeed(prop) * density(prop);
    m_rate_esti = A_th * sqrt(2 * k(prop) * pressure(prop) * density_0 * ((pressure_0 / pressure(prop)) ^ ((k(prop) - 1) / k(prop)) - 1) * (pressure(prop) / pressure_0) ^ (1 / k(prop)) / (k(prop) - 1));
    Err = (m_rate_calc - m_rate_esti)/m_rate_calc;

    %checking for roots
    if Err < 0;
        pressure_min = pressure_th;
    else
        pressure_max = pressure_th;
    end
    if abs(Err) < 10e-6
        break;
    end
end
m_rate = m_rate_calc;
density_th = density(prop);
pressure_th = (pressure_min + pressure_max) / 2;
set(prop, 'P', pressure_th, 'S', entropy_0);
temperature_th = temperature(prop);


%   Exit
for i = 1
    set(prop, 'P', pressure_th, 'S', entropy_0);
    A_exI = A_ex * i;
    mach_min = 1;
    mach_max = 100;

    while 1
        mach = (mach_min + mach_max) / 2;
        Arearatio_calc = ((1 + (k(prop) - 1) / 2 * (mach) ^ 2) / ((k(prop) + 1) / 2)) ^ ((k(prop) + 1) / (2 * (k(prop) - 1))) / (mach);
        Arearatio_actu = A_exI / A_th;
        Err = (Arearatio_calc - Arearatio_actu)/Arearatio_actu;

        %checking for roots
        if Err < 0;
            mach_min = mach;
            %disp("-")
        else
            mach_max = mach;
            %disp("+")
        end
        if abs(Err) < 10e-6
            break;
        end
    end
    j = 0;
    pressure_ex = pressure_0 * (1 + (k(prop) - 1) / 2 * (mach) ^ 2) ^ (- k(prop) / (k(prop) - 1));
    set(prop, 'P', pressure_ex, 'S', entropy_0);
    temperature_ex = temperature(prop);
    velocity_ex = mach * soundspeed(prop);
    pressure_min = pressure_ex / 2;
    pressure_max = pressure_ex + pressure_min;
    while 1
        pressure_ex = (pressure_min + pressure_max) / 2;
        set(prop, 'P', pressure_ex, 'S', entropy_0);
        temperature_ex = temperature(prop);
        velocity_ex = mach * soundspeed(prop);
        m_rate_esti = A_ex * velocity_ex * density(prop);
        Err = (m_rate_calc - m_rate_esti)/m_rate_calc
        if Err > 0;
            pressure_min = pressure_ex;
        else
            pressure_max = pressure_ex;
        end
        if abs(Err) < 10e-6
            break;
        end
        j = j + 1
        if j > 100
            break
        end
    end
    kRecord(i) = k(prop);
    Velocity_exRecord(i) = velocity_ex;
end
pressure_ex / oneatm
pressure_ex / pressure_th

%Functions:

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
k = 1.2;
%k = cp_mass(gas)/cv_mass(gas);
end