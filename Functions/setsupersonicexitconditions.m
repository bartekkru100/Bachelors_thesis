function s_supersonicExit = setsupersonicexitconditions(gas, expansionRatio)
import Gas.*

% This function calculates the exit state assuming a completely supersonic
% flow in the diverging part of the nozzle. Iterating over pressure to
% satisfy mass conservation.

    

%--------------------------------------------------------------------------

% Output
setarearatioisentropic(gas, expansionRatio, "super")
s_supersonicExit = State(gas);
end