function s_supersonicExit = setsupersonicexitconditions(gas, expansionRatio)
import Gas.*

s_throat = gas.sonic;

% This function calculates the exit state assuming a completely supersonic
% flow in the diverging part of the nozzle. Iterating over pressure to
% satisfy mass conservation.

setarearatioisentropic(gas, expansionRatio, "super");

%--------------------------------------------------------------------------

% Output

s_supersonicExit = State(gas);
end