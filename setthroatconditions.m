function gas = setthroatconditions(gas, s_stagnation)
import Gas.*

% This looks for sonic conditions (v = a)

% Using error correction method.

tolerance = 1e-12;
maxIterations = 100;
error_V_abs = -1;

numericalMethod = ErrorCorrectedMethod("setthroatconditions", tolerance, maxIterations);

while 1
    error_V_abs = error_V_abs + (1 - gas.Mach);

    pressure = s_stagnation.pressure * (1  + (s_stagnation.k - 1) / 2 * (1 + error_V_abs) ^ 2) ^ ((- s_stagnation.k) / (s_stagnation.k - 1));
    numericalMethod.findnewX(pressure);

    setpressureisentropic(gas, pressure, s_stagnation);
    error_V = (gas.soundspeed - gas.velocity) / gas.soundspeed; % Checking for convergence
    
    numericalMethod.updateXY(pressure, error_V);
    if numericalMethod.checkconvergence
        break;
    end
end

%--------------------------------------------------------------------------

% Output
s_throat = State(gas);
end