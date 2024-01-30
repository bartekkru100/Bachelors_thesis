function [pressure_shockAtExit] = findshockatexit_pressure(gas, s_maxVelocity, s_throat, s_stagnation, s_supersonicExit, expansionRatio)
import Gas.*

% Searching for an expansion ratio with shock at the exit. For a normal shock
% at the exit, the pressure after the shock will be equal to ambient
% pressure and velocities V_1 and V_2 (before and after the shock) fulfill
% the relation V_1 = V_throat ^ 2 / V_2. We look for such pairs of V_1 and
% V_2 such that mass is conserved across the shock. Expansion ratio is
% calculated from the ratio of mass flux at the throat to exit.

% Using Muller's method

tolerance = 1e-12;
iterationLimit = 100;

% Point 1
pressure(1) = s_maxVelocity.pressure;
error_M(1) = error_exit(gas, pressure(1), s_throat, s_supersonicExit);

% Point 2
pressure(2) = s_stagnation.pressure;
error_M(2) = error_exit(gas, pressure(2), s_throat, s_supersonicExit);

% Point 3
pressure(3) = s_throat.pressure;
error_M(3) = error_exit(gas, pressure(3), s_throat, s_supersonicExit);

numericalMethod = MullersMethod("checkfornormalshock_pressure", tolerance, iterationLimit, pressure, error_M, 'max');
numericalMethod.setX_min_max(s_maxVelocity.pressure, s_stagnation.pressure);

while 1
    pressure = numericalMethod.findnewX;

    error_M = error_exit(gas, pressure, s_throat, s_supersonicExit);

    numericalMethod.updateXY(pressure, error_M);
    if numericalMethod.checkconvergence
        break;
    end
end
%--------------------------------------------------------------------------

% Output

pressure_shockAtExit = pressure;

setstate(gas, s_supersonicExit);
end

%--------------------------------------------------------------------------

% Function to check mass conservation at given V_1

function error_M = error_exit(gas, pressure, s_throat, s_supersonic)
import Gas.*

setstate(gas, s_supersonic);
s_shock_1.massFlowFlux = gas.massFlowFlux;

s_shock_2.velocity = s_throat.velocity ^ 2 / gas.velocity;
s_shock_2.enthalpy = gas.totEnergy - s_shock_2.velocity ^ 2 / 2;
setstate(gas, 'H', s_shock_2.enthalpy, 'P', pressure, 'velocity', s_shock_2.velocity);

error_M = (s_shock_1.massFlowFlux - gas.massFlowFlux) / s_shock_1.massFlowFlux;
end