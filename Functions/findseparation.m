function [expansionRatio, s_separation]  = checkforseparation(gas, s_throat, s_atmo, separationTolerance)
import Gas.*

% Searching for separation using the criteria from the paper by
% Ralf H. Stark.

% Using Muller's method

tolerance = 1e-9;
iterationLimit = 100;

% Point 1
pressure(1) = s_atmo.pressure * 0.7;
error_Mach(1) = error_separation(gas, pressure(1), s_throat, s_atmo, separationTolerance);

% Point 2
pressure(2) = s_atmo.pressure * 0.45;
error_Mach(2) = error_separation(gas, pressure(2), s_throat, s_atmo, separationTolerance);

% Point 3
pressure(3) = s_atmo.pressure * 0.2;
error_Mach(3) = error_separation(gas, pressure(3), s_throat, s_atmo, separationTolerance);

numericalMethod = MullersMethod("checkforseparation", tolerance, 100, pressure, error_Mach, 'min');
numericalMethod.setX_min_max(0, s_atmo.pressure);

while 1
    pressure_0 = numericalMethod.findnewX();

    error_Mach_0 = error_separation(gas, pressure_0, s_throat, s_atmo, separationTolerance);

    numericalMethod.updateXY(pressure_0, error_Mach_0);
    if numericalMethod.checkconvergence;
        break;
    end
end

%--------------------------------------------------------------------------

% Output

expansionRatio = s_throat.massFlowFlux / gas.massFlowFlux;
s_separation = State(gas);
setstate(gas, s_throat);
end

%--------------------------------------------------------------------------

% Function that checks for the above mentioned criteria

function error_Mach = error_separation(gas, pressure, s_throat, s_atmo, separationTolerance)
import Gas.*
setpressureisentropic(gas, pressure, s_throat);

Mach_calc = s_atmo.pressure * (1 - separationTolerance) * pi / (3 * pressure);

error_Mach = (gas.Mach - Mach_calc) / gas.Mach;
end