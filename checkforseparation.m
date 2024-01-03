function [expansionRatio, distanceFromSeparation, s_separation]  = checkforseparation(gas, s_throat, s_exit, s_atmo)
import Gas.*

Tolerance = 1e-12;
% Using Muller's method
% Point 1
pressure_1 = s_atmo.pressure * 0.7;
error_Mach_1 = separationerror(gas, pressure_1, s_exit, s_atmo);

% Point 2
pressure_2 = s_atmo.pressure * 0.45;
error_Mach_2 = separationerror(gas, pressure_2, s_exit, s_atmo);

% Point 3
pressure_3 = s_atmo.pressure * 0.2;
error_Mach_3 = separationerror(gas, pressure_3, s_exit, s_atmo);

n = 0;
while 1
    n = n + 1;
    if n > 100
        disp("convergence failed in checkforseparation at error = " + error_Mach_1)
        break;
    end

    q = (pressure_1 - pressure_2) / (pressure_2 - pressure_3);
    a = q * error_Mach_1 - q * (1 + q) * error_Mach_2 + q ^ 2 * error_Mach_3;
    b = (2 * q + 1) * error_Mach_1 - (1 + q) ^ 2 * error_Mach_2 + q ^ 2 * error_Mach_3;
    c = (1 + q) * error_Mach_1;
    sqrtDelta = sqrt(b ^ 2 - 4 * a * c);
    denom(1) = (b + sqrtDelta);
    denom(2) = (b - sqrtDelta);
    pressure_0 = min(pressure_1 - (pressure_1 - pressure_2) * (2 * c) ./ denom);
    pressure_0 = real(pressure_0);
    if pressure_0 <= 0
        pressure_0 = s_atmo.pressure * 0.2 / n;
    end
    error_Mach_0 = separationerror(gas, pressure_0, s_exit, s_atmo);

    pressure_3 = pressure_2;
    error_Mach_3 = error_Mach_2;
    pressure_2 = pressure_1;
    error_Mach_2 = error_Mach_1;
    pressure_1 = pressure_0;
    error_Mach_1 = error_Mach_0;

    if abs(error_Mach_1) < Tolerance
        break;
    end
end
distanceFromSeparation = (s_exit.pressure - gas.pressure) / gas.pressure;
expansionRatio = s_throat.massFlowFlux / gas.massFlowFlux;
s_separation = State(gas);
s_separation.pressure / s_atmo.pressure
setstate(gas, s_exit);
end

function error_Mach = separationerror(gas, pressure, s_exit, s_atmo)
import Gas.*
setpressureisentropic(gas, pressure, s_exit);

Mach_calc = s_atmo.pressure * pi / (3 * pressure);

error_Mach = (gas.Mach - Mach_calc) / gas.Mach;
end