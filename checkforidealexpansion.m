function expansionRatio_ideal = checkforidealexpansion(gas, s_throat, s_atmo)
import Gas.*

% Searching for an expansion ratio with shock at the throat.

setpressureisentropic(gas, s_atmo.pressure);

expansionRatio_ideal = s_throat.massFlowFlux / gas.massFlowFlux;
setstate(gas, s_throat);
end