function expansionRatio_ideal = findidealexpansion(gas, s_throat, s_atmo)
import Gas.*

% Searching for an expansion ratio with shock at the throat or ideal expansion.

setpressureisentropic(gas, s_atmo.pressure);

expansionRatio_ideal = s_throat.massFlowFlux / gas.massFlowFlux;
setstate(gas, s_throat);
end