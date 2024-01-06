function [maxExpansion_shock, s_shockAtExit_1] = checkfornormalshock(gas, s_throat, s_atmo)
import State.*
import Gas.*


s_maxVelocity = gas.maxvelocity;
tolerance = 1e-4;

% searching for expansion ratio with shock at the exit

velocity(1) = s_maxVelocity.velocity;
error_M(1) = exiterror(gas, velocity(1), s_throat, s_atmo);
velocity(2) = s_throat.velocity;
error_M(2) = exiterror(gas, velocity(2), s_throat, s_atmo);

i = 1;
for vel = linspace(2232.52600268874, 2099.68972035407, 50)
    err(i, 1) = vel;
    err(i, 2) = exiterror(gas, vel, s_throat, s_atmo);
    err(i, 3) = s_throat.massFlowFlux / gas.massFlowFlux;
    i = i + 1
end

numericalMethod = BisectionMethod("checkfornormalshock_stage_1", tolerance, 100, velocity, error_M);
while 1
    velocity = numericalMethod.findnewX;
    error_M = exiterror(gas, velocity, s_throat, s_atmo);

    numericalMethod.updateXY(velocity, error_M);
    if numericalMethod.checkconvergence();
        break;
    end
end

tolerance = 1e-12;

velocity(1) = numericalMethod.get.X(1);
error_M(1) = numericalMethod.get.Y(1);

velocity(3) = numericalMethod.get.X(2);
error_M(3) = numericalMethod.get.Y(2);

velocity(2) = (velocity(1) + velocity(3)) / 2;
error_M(2) = exiterror(gas, velocity(2), s_throat, s_atmo);

numericalMethod = MullersMethod("checkfornormalshock_stage_2", tolerance, 100, velocity, error_M, 'max')

while 1
    velocity = numericalMethod.findnewX

    error_M = exiterror(gas, velocity, s_throat, s_atmo);

    numericalMethod.updateXY(velocity, error_M);
    if numericalMethod.checkconvergence
        break;
    end
end

s_shockAtExit_2 = State(gas);
setstate(gas, s_throat);
setvelocityisentropic(gas, velocity);
s_shockAtExit_1 = State(gas);

maxExpansion_shock = s_throat.massFlowFlux / gas.massFlowFlux;

% searching for expansion ratio with shock at the throat

setstate(gas, s_throat);
setpressureisentropic(gas, s_atmo.pressure, s_throat);
s_shockAtThroat = State(gas);

optimalExpansion_shock = s_throat.massFlowFlux / gas.massFlowFlux;
end













function error_M = exiterror(gas, velocity, s_throat, s_atmo)
import Gas.*

setstate(gas, s_throat);
setvelocityisentropic(gas, velocity);
s_shock_1 = State(gas);

s_shock_2.velocity = s_throat.velocity ^ 2 / velocity;
s_shock_2.enthalpy = s_shock_1.totEnergy - s_shock_2.velocity ^ 2 / 2;
setstate(gas, 'H', s_shock_2.enthalpy, 'P', s_atmo.pressure, 'velocity', s_shock_2.velocity);

error_M = (s_shock_1.massFlowFlux - gas.massFlowFlux) / s_shock_1.massFlowFlux;
end