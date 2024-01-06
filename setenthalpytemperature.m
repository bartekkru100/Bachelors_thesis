function setenthalpytemperature(gas, enthalpy_2, temperature_2)
import Gas.*

s_state_1 = State(gas);
Tolerance = 1e-12;

% Point 1
pressure(1) = 1;
error_H(1) = enthalpyerror(gas, enthalpy_2, pressure(1), temperature_2);

% Point 2
pressure(2) = s_state_1.pressure;
error_H(2) = enthalpyerror(gas, enthalpy_2, pressure(2), temperature_2);

% Point 3
pressure(3) = s_state_1.pressure * 2;
error_H(3) = enthalpyerror(gas, enthalpy_2, pressure(3), temperature_2);

n = 0;
while 1
    n = n + 1;
    if n > 100
        disp("convergence failed in setenthalpyisentropic with error = " + error_H(1));
        break;
    end

    pressure_0 = MullersMethod.findnewx(pressure, error_H, 'min');

    if pressure_0 <= 0
        pressure_0 = min([pressure(1), pressure(2), pressure(3)]) / 2;
    end
    error_H_0 = enthalpyerror(gas, enthalpy_2, pressure_0, state_1);

    [pressure, error_H] = MullersMethod.updatexy(pressure, error_H, pressure_0, error_H_0);

    if abs(error_H(1)) < Tolerance
        break;
    end
end

end

function error_H = enthalpyerror(gas, enthalpy_2, pressure, temperature_2)
set(gas.solution, 'T', temperature_2, 'P', pressure);
error_H = (enthalpy_2 - gas.enthalpy) / enthalpy_2;
end