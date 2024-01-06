function gas = setthroatconditions(gas)
import Gas.*

% Declaration of states.
state_1 = State(gas);
s_stagnation = gas.stagnation;

% Error declarations.
tolerance = 1e-12;
error_V = tolerance * 10;
error_V_abs = 0;

% Initial estimation of critical pressure.

numericalMethod = ErrorCorrectedMethod("setthroatconditions", tolerance, 100);

while 1
    error_V_abs = error_V_abs + (1 - gas.Mach);

    pressure = s_stagnation.pressure * (1  + (s_stagnation.k - 1) / 2 * (1 + error_V_abs) ^ 2) ^ ((- s_stagnation.k) / (s_stagnation.k - 1));
    numericalMethod.findnewX(pressure);

    setpressureisentropic(gas, pressure, state_1);
    error_V = (gas.soundspeed - gas.velocity) / gas.soundspeed; % Checking for convergence
    
    numericalMethod.updateXY(pressure, error_V);
    if numericalMethod.checkconvergence;
        break;
    end
end
s_throat = State(gas);
end