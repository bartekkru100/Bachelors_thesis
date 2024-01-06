function [expansionRatio, flowState, s_shock_1, s_shock_2] = shockposition(gas, s_atmo, s_supersonicExit, s_throat)
import State.*
import Gas.*

s_maxVelocity = gas.maxvelocity;
s_stagnation1 = gas.stagnation;
setstate(gas, 'P', s_atmo.pressure, 'H', s_stagnation1.enthalpy, 'velocity', 0);
setstate(gas, s_throat);
tolerance = 1e-12;

% Using Muller's method

% Point 1
velocity(1) = 0;
error_M(1) = exiterror(gas, velocity(1), s_throat, s_supersonicExit, s_atmo);

% Point 2
velocity(2) = s_throat.velocity / 2;
error_M(2) = exiterror(gas, velocity(2), s_throat, s_supersonicExit, s_atmo);

% Point 3
velocity(3) = s_throat.velocity;
error_M(3) = exiterror(gas, velocity(3), s_throat, s_supersonicExit, s_atmo);

numericalMethod = MullersMethod("shockposition", tolerance, 100, velocity, error_M, 'max');

while 1

    velocity_0 = numericalMethod.findnewX;

    if velocity_0 >= s_throat.velocity
        velocity_0 = s_throat.velocity * (numericalMethod.iteration + 1) / (numericalMethod.iteration + 2);
        if velocity(1) >= s_throat.velocity
            %error("failed to find shocked exit condition")
            break;
        end
    end
    error_M_0 = exiterror(gas, velocity_0, s_throat, s_supersonicExit, s_atmo);

    numericalMethod.updateXY(velocity_0, error_M_0);
    if numericalMethod.checkconvergence
        break;
    end
end

s_shockExit = State(gas);
% Using Muller's method

% Point 1
setstate(gas, s_throat);
setpressureisentropic(gas, s_atmo.pressure / s_maxVelocity.Mach);
velocity(1) = gas.velocity;
error_M(1) = shockerror(gas, velocity(1), s_throat, s_shockExit);

% Point 2
velocity(2) = (velocity(1) + s_maxVelocity.velocity) / 2;
error_M(2) = shockerror(gas, velocity(2), s_throat, s_shockExit);

% Point 3
velocity(3) = (velocity(1) + s_maxVelocity.velocity * 2) / 3;
error_M(3) = shockerror(gas, velocity(3), s_throat, s_shockExit);

velocity_fallback = velocity(3);

numericalMethod = MullersMethod("shockposition", tolerance, 100, velocity, error_M, 'max')

while 1

    velocity_0 = numericalMethod.findnewX;

    if velocity_0 >= s_maxVelocity.velocity
        velocity_0 = (velocity_fallback + s_maxVelocity.velocity * n) / (n + 1);
    end
    [error_M_0, s_shock_1, s_shock_2] = shockerror(gas, velocity_0, s_throat, s_shockExit);

    numericalMethod.updateXY(velocity_0, error_M_0);
    if numericalMethod.checkconvergence
        break;
    end
end
if velocity_0 > s_supersonicExit.velocity
    flowState = 'normal shock in the nozzle';
    setstate(gas, s_shockExit);
    expansionRatio = s_throat.massFlowFlux / gas.massFlowFlux;
else
    flowState = 'overexpanded';
    setstate(gas, s_supersonicExit);
    expansionRatio = s_throat.massFlowFlux / s_supersonicExit.massFlowFlux;
end
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