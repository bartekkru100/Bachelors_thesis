function [expansionRatio_shockAtExit, s_shockAtExit_1, s_shockAtExit_2, shockFound] = checkforshockatexit(gas, s_throat, s_atmo, s_maxVelocity, min_max)
import Gas.*

% Searching for an expansion ratio with shock at the exit. For a normal shock
% at the exit, the pressure after the shock will be equal to ambient
% pressure and velocities V_1 and V_2 (before and after the shock) fulfill
% the relation V_1 = V_throat ^ 2 / V_2. We look for such pairs of V_1 and
% V_2 such that mass is conserved across the shock. Expansion ratio is
% calculated from the ratio of mass flux at the throat to exit.

% Bisection is used to narrow down the search range quickly, then Muller's method
% is used to look for the root.

% Bisection method

tolerance = 1e-4;
iterationLimit = 10;

% Point 1
s_shockAtExit_1.pressure(1) = s_maxVelocity.pressure;
error_M(1) = exiterror(gas, s_shockAtExit_1.pressure(1), s_throat, s_atmo);

% Point 2
s_shockAtExit_1.pressure(2) = s_throat.pressure;
error_M(2) = exiterror(gas, s_shockAtExit_1.pressure(2), s_throat, s_atmo);

numericalMethod = BisectionMethod("checkfornormalshock_stage_1", tolerance, iterationLimit, s_shockAtExit_1.pressure, error_M);
numericalMethod.disablewarnings;

while 1
    s_shockAtExit_1.pressure = numericalMethod.findnewX;

    error_M = exiterror(gas, s_shockAtExit_1.pressure, s_throat, s_atmo);

    numericalMethod.updateXY(s_shockAtExit_1.pressure, error_M);
    if numericalMethod.checkconvergence
        break;
    end
end

%--------------------------------------------------------------------------

% Muller's method

tolerance = 1e-12;
iterationLimit = 100;

% Point 1
s_shockAtExit_1.pressure(1) = numericalMethod.get.X(1);
error_M(1) = numericalMethod.get.Y(1);

% Point 2
s_shockAtExit_1.pressure(3) = numericalMethod.get.X(2);
error_M(3) = numericalMethod.get.Y(2);

% Point 3
s_shockAtExit_1.pressure(2) = (s_shockAtExit_1.pressure(1) + s_shockAtExit_1.pressure(3)) / 2;
error_M(2) = exiterror(gas, s_shockAtExit_1.pressure(2), s_throat, s_atmo);

numericalMethod = MullersMethod("checkfornormalshock_stage_2", tolerance, iterationLimit, s_shockAtExit_1.pressure, error_M, min_max);
numericalMethod.setX_min_max(s_maxVelocity.pressure, s_throat.pressure);
numericalMethod.disablewarnings;

while 1
    s_shockAtExit_1.pressure = numericalMethod.findnewX;

    error_M = exiterror(gas, s_shockAtExit_1.pressure, s_throat, s_atmo);

    numericalMethod.updateXY(s_shockAtExit_1.pressure, error_M);
    if numericalMethod.checkconvergence
        if numericalMethod.hasfailed
            shockFound = false;
        else
            shockFound = true;
        end
        break;
    end
end

%--------------------------------------------------------------------------

% Output

s_shockAtExit_2 = State(gas);

setstate(gas, s_throat);
setpressureisentropic(gas, s_shockAtExit_1.pressure);
s_shockAtExit_1 = State(gas);

expansionRatio_shockAtExit = s_throat.massFlowFlux / gas.massFlowFlux;

setstate(gas, s_throat);
end

%--------------------------------------------------------------------------

% Function to check mass conservation at given V_1

function error_M = exiterror(gas, pressure, s_throat, s_atmo)
import Gas.*

setstate(gas, s_throat);
setpressureisentropic(gas, pressure, s_throat);
s_shock_1.massFlowFlux = gas.massFlowFlux;

s_shock_2.velocity = s_throat.velocity ^ 2 / gas.velocity;
s_shock_2.enthalpy = gas.totEnergy - s_shock_2.velocity ^ 2 / 2;
setstate(gas, 'H', s_shock_2.enthalpy, 'P', s_atmo.pressure, 'velocity', s_shock_2.velocity);

error_M = (s_shock_1.massFlowFlux - gas.massFlowFlux) / s_shock_1.massFlowFlux;
end