% Author: Bartosz Kruszona, 2023

% Injection
% Cantera sets the injection composition, pressure and temperature
function injection(prop, state)
import State.*
setgasstate(prop, state);
end