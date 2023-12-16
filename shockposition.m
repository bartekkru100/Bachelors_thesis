function [expansionRatio, flowState] = shockposition(gas, s_atmo, s_exit, s_throat)
import State.*
import Gas.*

s_stagnation1 = gas.stagnation;
setstate(gas, 'P', s_atmo.pressure, 'H', s_stagnation1.enthalpy, 'velocity', 0);
s_stagnation2 = State(gas);
setstate(gas, s_throat);
Tolerance = 1e-12;

tic
% Using Muller's method

% Point 1
velocity_1 = s_throat.velocity;
error_M_1 = setshockstate(gas, velocity_1, s_throat, s_stagnation2);

% Point 2
velocity_2 = s_throat.velocity * 2;
error_M_2 = setshockstate(gas, velocity_2, s_throat, s_stagnation2);

% Point 3
velocity_3 = s_exit.velocity;
error_M_3 = setshockstate(gas, velocity_3, s_throat, s_stagnation2);

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

    if velocity_0 >= s_exit.velocity
        velocity_0 = max([velocity_1, velocity_2, velocity_3, s_exit.velocity]);
        if velocity_1 >= s_exit.velocity
            if s_exit.pressure > s_atmo.pressure
                flowState = 'underexpanded'
            else
                flowState = 'overexpanded'
            end
            break;
        end
    end
    error_M_0 = setshockstate(gas, velocity_0, s_throat, s_stagnation2);

    velocity_3 = velocity_2;
    error_M_3 = error_M_2;
    velocity_2 = velocity_1;
    error_M_2 = error_M_1;
    velocity_1 = velocity_0;
    error_M_1 = error_M_0;

    if abs(error_M_1) < Tolerance
        flowState = 'normal shock in the nozzle'
        break;
    end
end
toc


%{
tic
% Using Muller's method

% Point 1
velocity_1 = s_throat.velocity / 2;
error_M_1 = setshockstate(gas, velocity_1, s_throat, s_stagnation2);

% Point 2
velocity_2 = s_exit.velocity;
error_M_2 = setshockstate(gas, velocity_2, s_throat, s_stagnation2);

n = 0;
while 1
    n = n + 1;
    if n > 100
        disp("convergence failed")
        break;
    end

    velocity_0 = (velocity_2 * error_M_1 - velocity_1 * error_M_2 ) / (error_M_1 - error_M_2)

    error_M_0 = setshockstate(gas, velocity_0, s_throat, s_stagnation2);

    velocity_2 = velocity_1;
    error_M_2 = error_M_1;
    velocity_1 = velocity_0;
    error_M_1 = error_M_0;

    if abs(error_M_1) < Tolerance
        break;
    end
end
toc
%}
expansionRatio = s_throat.massFlowFlux / gas.massFlowFlux;
end

function error_M = setshockstate(gas, velocity, s_throat, s_stagnation2)
import Gas.*
setstate(gas, s_throat);
setvelocityisentropic(gas, velocity);
s_shock1 = State(gas);

setstate(gas, s_stagnation2);
s_shock2.velocity = s_throat.velocity ^ 2 / s_shock1.velocity;
setvelocityisentropic(gas, s_shock2.velocity);
s_shock2 = State(gas);

error_M = (s_shock1.massFlowFlux - s_shock2.massFlowFlux) / s_shock1.massFlowFlux;
end