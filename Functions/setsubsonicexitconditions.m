function [s_subsonicExit, s_throat, s_chamber] = setsubsonicexitconditions(gas, s_atmo, expansionRatio, contractionRatio)
import Gas.*

% Finds the chamber and exit conditions for a subsonic flow

s_stagnation = gas.stagnation(false);
setpressureisentropic(gas, s_atmo.pressure, s_stagnation);
s_subsonicExit = State(gas);
setarearatioisentropic(gas, s_subsonicExit.areaRatio / expansionRatio, "sub");
s_throat = State(gas);
setarearatioisentropic(gas, s_throat.areaRatio * contractionRatio, "sub");
s_chamber = State(gas);

end