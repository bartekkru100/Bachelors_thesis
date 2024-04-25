function pressure_separation  = findseparation_pressure(gas, s_supersonicExit, separationTolerance)
import Gas.*
pressure_separation = 3 * s_supersonicExit.pressure * s_supersonicExit.Mach / (pi * (1 - separationTolerance));
end