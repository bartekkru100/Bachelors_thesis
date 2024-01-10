function [expansionRatio_shockAtExit, s_shockAtExit_1, s_shockAtExit_2, flowState] = checkfornormalshock(gas, s_throat, s_atmo, s_maxVelocity)
import Gas.*

% Searching for expansion ratio with shock at the exit. For a normal shock
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
s_shockAtExit_1.velocity(1) = s_maxVelocity.velocity;
error_M(1) = exiterror(gas, s_shockAtExit_1.velocity(1), s_throat, s_atmo);

% Point 2
s_shockAtExit_1.velocity(2) = s_throat.velocity;
error_M(2) = exiterror(gas, s_shockAtExit_1.velocity(2), s_throat, s_atmo);

numericalMethod = BisectionMethod("checkfornormalshock_stage_1", tolerance, iterationLimit, s_shockAtExit_1.velocity, error_M);
numericalMethod.disablewarnings;

    tic
while 1
    s_shockAtExit_1.velocity = numericalMethod.findnewX;

    error_M = exiterror(gas, s_shockAtExit_1.velocity, s_throat, s_atmo);

    numericalMethod.updateXY(s_shockAtExit_1.velocity, error_M);
    if numericalMethod.checkconvergence
        break;
    end
end
    toc

%--------------------------------------------------------------------------

% Muller's method

tolerance = 1e-12;
iterationLimit = 100;

% Point 1
s_shockAtExit_1.velocity(1) = numericalMethod.get.X(1);
error_M(1) = numericalMethod.get.Y(1);

% Point 2
s_shockAtExit_1.velocity(3) = numericalMethod.get.X(2);
error_M(3) = numericalMethod.get.Y(2);

% Point 3
s_shockAtExit_1.velocity(2) = (s_shockAtExit_1.velocity(1) + s_shockAtExit_1.velocity(3)) / 2;
error_M(2) = exiterror(gas, s_shockAtExit_1.velocity(2), s_throat, s_atmo);

numericalMethod = MullersMethod("checkfornormalshock_stage_2", tolerance, iterationLimit, s_shockAtExit_1.velocity, error_M, 'max');
numericalMethod.setX_min_max(s_throat.velocity, s_maxVelocity.velocity);
numericalMethod.disablewarnings;

while 1
    s_shockAtExit_1.velocity = numericalMethod.findnewX;

    error_M = exiterror(gas, s_shockAtExit_1.velocity, s_throat, s_atmo);

    numericalMethod.updateXY(s_shockAtExit_1.velocity, error_M);
    if numericalMethod.checkconvergence
        if numericalMethod.hasfailed
            flowState = "No exit shocks possible";
        else
            flowState = "Shock at the exit found";
        end
        break;
    end
end

%--------------------------------------------------------------------------

% Output

s_shockAtExit_2 = State(gas);

setstate(gas, s_throat);
setvelocityisentropic(gas, s_shockAtExit_1.velocity);
s_shockAtExit_1 = State(gas);

expansionRatio_shockAtExit = s_throat.massFlowFlux / gas.massFlowFlux;

setstate(gas, s_throat);
end

%--------------------------------------------------------------------------

% Function to check mass conservation at given V_1

function error_M = exiterror(gas, velocity, s_throat, s_atmo)
import Gas.*

setstate(gas, s_throat);
setvelocityisentropic(gas, velocity);
s_shock_1.massFlowFlux = gas.massFlowFlux;

s_shock_2.velocity = s_throat.velocity ^ 2 / velocity;
s_shock_2.enthalpy = gas.totEnergy - s_shock_2.velocity ^ 2 / 2;
setstate(gas, 'H', s_shock_2.enthalpy, 'P', s_atmo.pressure, 'velocity', s_shock_2.velocity);

error_M = (s_shock_1.massFlowFlux - gas.massFlowFlux) / s_shock_1.massFlowFlux;
end