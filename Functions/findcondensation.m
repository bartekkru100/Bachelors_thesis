function [expansionRatio_condensation, condensationFound] = findcondensation(gas)
import Gas.*

% Bisection method:

tolerance = 1e-9;
iterationLimit = 100;

pressure(1) = gas.s_stagnation.pressure;
setpressureisentropic(gas, pressure(1));
error_P(1) = 1;
if gas.hasCondensation
    error_P(1) = - error_P(1);
end

pressure(2) = gas.maxvelocity.pressure;
setpressureisentropic(gas, pressure(2));
error_P(2) = 1;
if gas.hasCondensation
    error_P(2) = - error_P(2);
end

if error_P(1) ~= error_P(2)
    numericalMethod = BisectionMethod("findcondensation", tolerance, iterationLimit, pressure, error_P);

    while 1
        pressure = numericalMethod.findnewX;

        error_P = error_condensation(gas, numericalMethod);
        numericalMethod.updateXY(pressure, error_P);
        if numericalMethod.checkconvergence
            break;
        end
    end

    condensationFound = true;

else
    condensationFound = false;

end

%--------------------------------------------------------------------------

% Output

expansionRatio_condensation = gas.s_sonic.massFlowFlux / gas.massFlowFlux;
setstate(gas, gas.s_sonic);
end

%--------------------------------------------------------------------------

% This function calculates the pressure error for the condensation point.
% The error is defined as the distance between bisection points, if
% condensation happens at the pressure midpoint, the sign is changed to negative.

function error_P = error_condensation(gas, numericalMethod)
import Gas.*

error_P = abs(numericalMethod.get.X(1) - numericalMethod.get.X(2));
setpressureisentropic(gas, mean(numericalMethod.get.X));
if gas.hasCondensation
    error_P = - error_P;
end
end