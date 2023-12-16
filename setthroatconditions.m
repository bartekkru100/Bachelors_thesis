function gas = setthroatconditions(gas)
import Gas.*

% Declaration of states.
s_stagnation = gas.stagnation;
s_throat = State(gas);

% Error declarations.
Tolerance = 1e-12;
error_V = Tolerance * 10;
error_V_old = error_V;
error_V_abs = 0;

% Initial estimation of critical pressure.

n = 0;
while abs(error_V) > Tolerance
    if abs(error_V_old / error_V) < 1
        n = n + 1;
        if n > 2
            break;
        end
    else
        n = 0;
    end
    error_V_old = error_V;

    error_V_abs = error_V_abs + (1 - s_throat.Mach);
    s_throat.pressure = s_stagnation.pressure * (1  + (s_stagnation.k - 1) / 2 * (1 + error_V_abs) ^ 2) ^ ((- s_stagnation.k) / (s_stagnation.k - 1));
    setPressureIsentropic(gas, s_throat.pressure);
    s_throat = State(gas);

    error_V = (s_throat.soundspeed - s_throat.velocity) / s_throat.soundspeed; % Checking for convergence
end
end