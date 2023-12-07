% Author: Bartosz Kruszona
import State.*
import Gas.*

clc, clear

% Constraints
geometry_defined = false;
expansion_ratio = 80;
contraction_ratio = 5;
heat_power_area = 0;
thrust = 1;

% Atmospheric conditions
atmo = Gas('Air.yaml');
temperature_at = 293;
pressure_at = 0.1;


% Initial Injection
pressure_in = 300e5;
temperature_in = 293;
prop = Gas('R1highT.yaml');
%prop = Gas('GRI30.yaml');  
%set(prop, 'Y', 'H2:1,O2:6.03');
setState(prop, 'Y', 'CH4:1,O2:3.6');
%set(prop, 'X', 'N2:1');
%set(prop, 'Y', 'POSF7688:1,O2:2.72');
setState(prop, 'P', pressure_in, 'T', temperature_in);
s_injection = State(prop);
injection(prop, s_injection);

