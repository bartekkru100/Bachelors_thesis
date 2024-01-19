function [s_subsonicExit, s_chamber] = setsubsonicexitconditions(gas, s_chamber, s_sonic, s_atmo, expansionRatio, contractionRatio)
import Gas.*

% Finds the chamber and exit conditions for a subsonic flow

tolerance = 1e-12;
maxIterations = 100;

velocity(1) = 0;
error_M = error_exitchamber(gas, s_chamber, s_atmo, velocity, expansionRatio, contractionRatio);

velocity(2) = s_sonic.velocity / 2;
error_M = error_exitchamber(gas, s_chamber, s_atmo, velocity, expansionRatio, contractionRatio);

velocity(3) = s_sonic.velocity;
error_M = error_exitchamber(gas, s_chamber, s_atmo, velocity, expansionRatio, contractionRatio);

numericalMethod = MullersMethod("setsubsonicexitconditions", tolerance, maxIterations, velocity, error_M, 'min');
numericalMethod.setX_min_max(0, s_sonic.velocity);

while 1
    velocity = numericalMethod.findnewX;

    error_M = error_exitchamber(gas, s_chamber, s_atmo, velocity, expansionRatio, contractionRatio);

    numericalMethod.updateXY(velocity, error_M);
    if numericalMethod.checkconvergence
        break;
    end
end

%--------------------------------------------------------------------------

% Output

if nargout >= 1
    s_subsonicExit = State(gas);
end
if nargout == 2
    setpressureisentropic(gas, s_chamber.pressure, s_subsonicExit);
    s_chamber = State(gas);
end

end

function error_M = error_exitchamber(gas, s_chamber, s_atmo, velocity, expansionRatio, contractionRatio)
import Gas.*

s_chamber.velocity = velocity;
setstate(gas, s_chamber);
s_chamber.massFlowFlux = gas.massFlowFlux;

setpressureisentropic(gas, s_atmo.pressure, s_chamber);
s_subsonicExit.massFlowFlux = gas.massFlowFlux;

error_M = (s_subsonicExit.massFlowFlux - s_chamber.massFlowFlux * contractionRatio / expansionRatio) / s_subsonicExit.massFlowFlux;
end