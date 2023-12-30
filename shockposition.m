function [expansionRatio, flowState, s_shock_1, s_shock_2] = shockposition(gas, s_atmo, s_supersonicExit, s_throat)
import State.*
import Gas.*

if s_supersonicExit.pressure > s_atmo.pressure
    flowState = 'underexpanded';
    s_shock_1 = 0;
    s_shock_2 = 0;
    expansionRatio = s_throat.massFlowFlux / s_supersonicExit.massFlowFlux;
    return;
end

s_stagnation1 = gas.stagnation;
setstate(gas, 'P', s_atmo.pressure, 'H', s_stagnation1.enthalpy, 'velocity', 0);
setstate(gas, s_throat);
Tolerance = 1e-12;

% Using Muller's method

% Point 1
velocity_1 = 0;
error_M_1 = exiterror(gas, velocity_1, s_throat, s_supersonicExit, s_atmo);

% Point 2
velocity_2 = s_throat.velocity / 2;
error_M_2 = exiterror(gas, velocity_2, s_throat, s_supersonicExit, s_atmo);

% Point 3
velocity_3 = s_throat.velocity;
error_M_3 = exiterror(gas, velocity_3, s_throat, s_supersonicExit, s_atmo);

n = 0;
while 1
    n = n + 1;
    if n > 100
        disp("convergence failed")
        break;
    end

    q = (velocity_1 - velocity_2) / (velocity_2 - velocity_3);
    a = q * error_M_1 - q * (1 + q) * error_M_2 + q ^ 2 * error_M_3;
    b = (2 * q + 1) * error_M_1 - (1 + q) ^ 2 * error_M_2 + q ^ 2 * error_M_3;
    c = (1 + q) * error_M_1;
    sqrtDelta = sqrt(b ^ 2 - 4 * a * c);
    denom(1) = (b + sqrtDelta);
    denom(2) = (b - sqrtDelta);
    velocity_0 = max(velocity_1 - (velocity_1 - velocity_2) * (2 * c) ./ denom);

    if velocity_0 >= s_throat.velocity
        velocity_0 = max([velocity_1, velocity_2, velocity_3, s_throat.velocity]);
        if velocity_1 >= s_throat.velocity
            break;
        end
    end
    error_M_0 = exiterror(gas, velocity_0, s_throat, s_supersonicExit, s_atmo);

    velocity_3 = velocity_2;
    error_M_3 = error_M_2;
    velocity_2 = velocity_1;
    error_M_2 = error_M_1;
    velocity_1 = velocity_0;
    error_M_1 = error_M_0;

    if abs(error_M_1) < Tolerance
        break;
    end
end

s_shockExit = State(gas);
% Using Muller's method

% Point 1
velocity_1 = s_throat.velocity;
error_M_1 = shockerror(gas, velocity_1, s_throat, s_shockExit);

% Point 2
velocity_2 = s_throat.velocity * 2;
error_M_2 = shockerror(gas, velocity_2, s_throat, s_shockExit);

% Point 3
velocity_3 = s_supersonicExit.velocity;
error_M_3 = shockerror(gas, velocity_3, s_throat, s_shockExit);

n = 0;
while 1
    n = n + 1;
    if n > 100
        disp("convergence failed")
        break;
    end

    q = (velocity_1 - velocity_2) / (velocity_2 - velocity_3);
    a = q * error_M_1 - q * (1 + q) * error_M_2 + q ^ 2 * error_M_3;
    b = (2 * q + 1) * error_M_1 - (1 + q) ^ 2 * error_M_2 + q ^ 2 * error_M_3;
    c = (1 + q) * error_M_1;
    sqrtDelta = sqrt(b ^ 2 - 4 * a * c);
    denom(1) = (b + sqrtDelta);
    denom(2) = (b - sqrtDelta);
    velocity_0 = max(velocity_1 - (velocity_1 - velocity_2) * (2 * c) ./ denom);

    if velocity_0 >= s_supersonicExit.velocity
        velocity_0 = max([velocity_1, velocity_2, velocity_3, s_supersonicExit.velocity]);
        if velocity_1 >= s_supersonicExit.velocity
            break;
        end
    end
    [error_M_0, s_shock_1, s_shock_2] = shockerror(gas, velocity_0, s_throat, s_shockExit);

    velocity_3 = velocity_2;
    error_M_3 = error_M_2;
    velocity_2 = velocity_1;
    error_M_2 = error_M_1;
    velocity_1 = velocity_0;
    error_M_1 = error_M_0;

    if abs(error_M_1) < Tolerance
        flowState = 'normal shock in the nozzle';
        expansionRatio = s_throat.massFlowFlux / gas.massFlowFlux;
        setstate(gas, s_shockExit);
        return;
    end
end
flowState = 'overexpanded';
setstate(gas, s_supersonicExit);
expansionRatio = s_throat.massFlowFlux / s_supersonicExit.massFlowFlux;
end


function error_M = exiterror(gas, velocity, s_throat, s_exit, s_atmo)
import Gas.*

setstate(gas, s_throat);
s_shockExit.enthalpy = s_throat.totEnergy - velocity ^ 2 / 2;
setstate(gas, 'H', s_shockExit.enthalpy, 'P', s_atmo.pressure, 'velocity', velocity);
s_shockExit = State(gas);

error_M = (s_exit.massFlowFlux - s_shockExit.massFlowFlux) / s_exit.massFlowFlux;
end

function [error_M, s_shock_1, s_shock_2] = shockerror(gas, velocity, s_throat, s_shockExit)
import Gas.*

setstate(gas, s_throat);
setvelocityisentropic(gas, velocity);
s_shock_1 = State(gas);

setstate(gas, s_shockExit);
s_shock_2.velocity = s_throat.velocity ^ 2 / s_shock_1.velocity;
setvelocityisentropic(gas, s_shock_2.velocity);
s_shock_2 = State(gas);

error_M = (s_shock_1.massFlowFlux - s_shock_2.massFlowFlux) / s_shock_1.massFlowFlux;
end