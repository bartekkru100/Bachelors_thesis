function [shockPosition, s_shock_1, s_shock_2, s_subsonicExit] = shockposition(gas, s_atmo, s_supersonicExit, s_throat)
import State.*
import Gas.*
i = 0;

setstate(gas, s_supersonicExit);

% Looking for subsonic velocity, for which the mass conservation will be
% maintained between sonic and subsonic flow at the exit. At exit, the
% pressure of subsonic flow is equal to the ambient pressure. Enthalpy is
% chosen so that total flow energy is conserved.

% Muller's method:

tolerance = 1e-9;
iterationLimit = 100;

% Point 1
s_subsonicExit.velocity(1) = 0;
error_M(1) = error_exit(gas, s_subsonicExit.velocity(1), s_supersonicExit, s_atmo);

% Point 2
s_subsonicExit.velocity(2) = s_throat.velocity / 2;
error_M(2) = error_exit(gas, s_subsonicExit.velocity(2), s_supersonicExit, s_atmo);

% Point 2
s_subsonicExit.velocity(3) = s_throat.velocity;
error_M(3) = error_exit(gas, s_subsonicExit.velocity(3), s_supersonicExit, s_atmo);

numericalMethod = MullersMethod("shockposition_stage_1", tolerance, iterationLimit, s_subsonicExit.velocity, error_M, 'max');

while 1
    s_subsonicExit.velocity = numericalMethod.findnewX;
    if s_subsonicExit.velocity < 0
        s_subsonicExit.velocity = 0;
    elseif s_subsonicExit.velocity > s_throat.velocity
        s_subsonicExit.velocity = s_throat.velocity;
    elseif isnan(s_subsonicExit.velocity)
        %error("Couldn't find a valid subsonic exit");
        break;
    end

    error_M = error_exit(gas, s_subsonicExit.velocity, s_supersonicExit, s_atmo);

    numericalMethod.updateXY(s_subsonicExit.velocity, error_M);
    if numericalMethod.checkconvergence
        break;
    end
end

error_M_test = error_M;

s_subsonicExit = State(gas);

%--------------------------------------------------------------------------

% Searching for pairs of V_1, V_2 (velocities before and after the shock)
% where V_1 = V_throat ^ 2 / V_2, for which mass is conserving assuming
% isentropic processes frm throat to the shock and then from the shock to
% the exit.

% Using bisection method initially, then switching to Muller's method.
% Using supersonic velocity as input and calculating subsonic velocity from
% it

% Bisection method:

tolerance = 1e-4;
iterationLimit = 10;

s_shock_1.velocity(1) = s_throat.velocity;
error_M(1) = error_shock(gas, s_shock_1.velocity(1), s_throat, s_subsonicExit);

s_shock_1.velocity(2) = s_supersonicExit.velocity;
error_M(2) = error_shock(gas, s_shock_1.velocity(2), s_throat, s_subsonicExit);

numericalMethod = BisectionMethod("shockposition_stage_2.1", tolerance, iterationLimit, s_shock_1.velocity, error_M);
numericalMethod.disablewarnings;

while 1
    s_shock_1.velocity = numericalMethod.findnewX;

    error_M = error_shock(gas, s_shock_1.velocity, s_throat, s_subsonicExit);
    numericalMethod.updateXY(s_shock_1.velocity, error_M);
    if numericalMethod.checkconvergence
        break;
    end
end

% Muller's method:

tolerance = 1e-9;
iterationLimit = 100;

% Point 1
s_shock_1.velocity(1) = numericalMethod.get.X(1);
error_M(1) = numericalMethod.get.Y(1);

% Point 2
s_shock_1.velocity(2) = mean(numericalMethod.get.X);
error_M(2) = error_shock(gas, s_shock_1.velocity(2), s_throat, s_subsonicExit);

% Point 2
s_shock_1.velocity(3) = numericalMethod.get.X(2);;
error_M(3) = numericalMethod.get.Y(2);

numericalMethod = MullersMethod("shockposition_stage_2.2", tolerance, iterationLimit, s_shock_1.velocity, error_M, 'max');

while 1
    s_shock_1.velocity = numericalMethod.findnewX;
    if isnan(s_shock_1.velocity)
        s_shock_1.velocity = max(numericalMethod.get.X) * (1 + (numericalMethod.get.Y_0_old));
        numericalMethod.setnewX(s_shock_1.velocity);
    end
    if s_shock_1.velocity > s_supersonicExit.velocity
        s_shock_1.velocity = s_supersonicExit.velocity;
    elseif s_shock_1.velocity < s_throat.velocity
        s_shock_1.velocity = s_throat.velocity;
    end

    error_M = error_shock(gas, s_shock_1.velocity, s_throat, s_subsonicExit);
    numericalMethod.updateXY(s_shock_1.velocity, error_M);
    if numericalMethod.checkconvergence
        break;
    end
end

%--------------------------------------------------------------------------

% Output

s_shock_2 = State(gas);

setstate(gas, s_throat);
setvelocityisentropic(gas, s_shock_1.velocity);
s_shock_1 = State(gas);

shockPosition = s_throat.massFlowFlux / gas.massFlowFlux;

setstate(gas, s_subsonicExit);
end

%--------------------------------------------------------------------------

% Function that find the subsonic gas state at the exit and checks for mass conservation.

function error_M = error_exit(gas, velocity, s_supersonicExit, s_atmo)
import Gas.*

setstate(gas, s_supersonicExit);
enthalpy = s_supersonicExit.totEnergy - velocity ^ 2 / 2;
setstate(gas, 'H', enthalpy, 'P', s_atmo.pressure, 'velocity', velocity);

error_M = (s_supersonicExit.massFlowFlux - gas.massFlowFlux) / s_supersonicExit.massFlowFlux;
end

%--------------------------------------------------------------------------

% Function that finds the state of the gas after the shock and checks for
% mass conservation.

function error_M = error_shock(gas, velocity, s_throat, s_subsonicExit)
import Gas.*

setstate(gas, s_throat);
setvelocityisentropic(gas, velocity);
s_shock_1.massFlowFlux = gas.massFlowFlux;

setstate(gas, s_subsonicExit);
s_shock_2.velocity = s_throat.velocity ^ 2 / velocity;
setvelocityisentropic(gas, s_shock_2.velocity);

error_M = (s_shock_1.massFlowFlux - gas.massFlowFlux) / s_shock_1.massFlowFlux;
end