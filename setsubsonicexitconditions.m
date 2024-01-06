function gas = setsubsonicexitconditions(gas, s_chamber, s_atmo)
import Gas.*

setstate(gas, s_chamber);
setpressureisentropic(gas, s_atmo.pressure, s_chamber);
end