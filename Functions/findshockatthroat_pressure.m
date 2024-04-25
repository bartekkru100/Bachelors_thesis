function pressure_shockAtThroat = findidealexpansion_pressure(gas, expansionRatio)
import Gas.*

% Returns an ambient pressure with shock at throat  for a given expansion ratio.

pressure_shockAtThroat = setarearatioisentropic(gas, expansionRatio, "sub").pressure;

end