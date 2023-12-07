function expansionRatio = shockposition(gas, s_stagnation, s_atmo, s_exit, s_throat)
import State.*
prop = gas;
s_shock1 = State(s_throat)
s_shock2 = State(s_throat)

velocity_min = s_throat.velocity;
velocity_max = s_exit.velocity;
for i = 1 : 30
    s_throat.velocity = (velocity_min + velocity_max) / 2

    mass_rate_1 = 
    mass_rate_2
    error_M;
    if abs(error_M) < 1e-6  % Checking for convergence
        break;
    end
    if error_M < 0 % Range redefinition
        velocity_min = s_throat.velocity;
    else
        velocity_max = s_throat.velocity;
    end
end

end