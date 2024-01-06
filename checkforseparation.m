function [expansionRatio, s_separation]  = checkforseparation(gas, s_throat, s_atmo, separationTolerance)
import Gas.*

tolerance = 1e-12;
% Using Muller's method
% Point 1
pressure(1) = s_atmo.pressure * 0.7;
error_Mach(1) = separationerror(gas, pressure(1), s_throat, s_atmo, separationTolerance);

% Point 2
pressure(2) = s_atmo.pressure * 0.45;
error_Mach(2) = separationerror(gas, pressure(2), s_throat, s_atmo, separationTolerance);

% Point 3
pressure(3) = s_atmo.pressure * 0.2;
error_Mach(3) = separationerror(gas, pressure(3), s_throat, s_atmo, separationTolerance);

numericalMethod = MullersMethod("checkforseparation", tolerance, 100, pressure, error_Mach, 'min');

n = 0;
while 1

    pressure_0 = numericalMethod.findnewX();

    if pressure_0 <= 0
        pressure_0 = s_atmo.pressure * 0.2 / n;
    end
    error_Mach_0 = separationerror(gas, pressure_0, s_throat, s_atmo, separationTolerance);

    numericalMethod.updateXY(pressure_0, error_Mach_0);
    if numericalMethod.checkconvergence;
        break;
    end
end
%distanceFromSeparation = (s_throat.pressure - gas.pressure) / gas.pressure;
expansionRatio = s_throat.massFlowFlux / gas.massFlowFlux;
s_separation = State(gas);
s_separation.pressure / s_atmo.pressure
setstate(gas, s_throat);
end

function error_Mach = separationerror(gas, pressure, s_throat, s_atmo, separationTolerance)
import Gas.*
setpressureisentropic(gas, pressure, s_throat);

Mach_calc = s_atmo.pressure * (1 - separationTolerance) * pi / (3 * pressure);

error_Mach = (gas.Mach - Mach_calc) / gas.Mach;
end