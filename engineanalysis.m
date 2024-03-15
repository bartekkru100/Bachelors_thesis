function solverOutputs = engineanalysis(solverInputs)

import Gas.*
import State.*

gas = solverInputs.gas;
atmo = solverInputs.atmo;
expansionRatio = solverInputs.expansionRatio;
contractionRatio = solverInputs.contractionRatio;
heatPowerArea = solverInputs.heatPowerArea;
thrust = solverInputs.thrust;

s_injection = State(gas);
setchamberconditions(gas, s_injection)
s_chamber = State(gas);

Tolerance = 1e-9;
massFlowFlux_old = 0;
for i = 1 : 20

    [s_throat, s_stagnation] = gas.sonic;

    error_M = (massFlowFlux_old - s_throat.massFlowFlux) / massFlowFlux_old;
    if abs(error_M) < Tolerance
        break;
    end

    setinjectionconditions(gas, s_injection, contractionRatio);
    s_injection = State(gas);

    setchamberconditions(gas, s_injection, heatPowerArea);
    s_chamber = State(gas);

    massFlowFlux_old = s_throat.massFlowFlux;
end
s_maxVelocity = gas.maxvelocity();

displaygasinfo(s_injection, "Injection", contractionRatio);
displaygasinfo(s_chamber, "Chamber", contractionRatio);
displaygasinfo(s_throat, "Throat", 1);

s_supersonicExit = setsupersonicexitconditions(gas, s_throat, expansionRatio);
%findshockatexit_pressure(gas, s_maxVelocity, s_throat, s_stagnation, s_supersonicExit, expansionRatio)

expansionRatio_condensation = findcondensation(gas, s_throat, s_stagnation, s_maxVelocity);

if s_stagnation.pressure > atmo.pressure
    [expansionRatio_shockAtExit, exitShockfound] = findshockatexit(gas, s_throat, atmo, s_maxVelocity, 'min');
    expansionRatio_shockAtThroat = findidealexpansion(gas, s_throat, atmo);
    s_supersonicExit = setsupersonicexitconditions(gas, s_throat, expansionRatio);

    if s_throat.pressure < atmo.pressure

        if expansionRatio < expansionRatio_shockAtThroat
            [s_exit, s_chamber] = setsubsonicexitconditions(gas, s_chamber, s_throat, atmo, expansionRatio, contractionRatio);
            flowState = "Fully subsonic";

        elseif expansionRatio <= expansionRatio_shockAtExit
            [expansionRatio_shockAtExit_2, exitShockfound] = findshockatexit(gas, s_throat, atmo, s_maxVelocity, max);

            if expansionRatio >= expansionRatio_shockAtExit_2
                flowState = "Fully supersonic";

            else
                [shockPosition, s_shock_1, s_shock_2] = shockposition(gas, atmo, s_supersonicExit, s_throat);
                flowState = "Normal shock inside the nozzle";

            end
        else
            [shockPosition, s_shock_1, s_shock_2] = shockposition(gas, atmo, s_supersonicExit, s_throat);
            flowState = "Normal shock inside the nozzle";

        end
    else
        if expansionRatio <= expansionRatio_shockAtExit
            flowState = "Fully supersonic";

        else
            [shockPosition, s_shock_1, s_shock_2] = shockposition(gas, atmo, s_supersonicExit, s_throat);
            flowState = "Normal shock inside the nozzle";

        end
    end
else
    flowState = "No flow possible";
end
s_exit = State(gas);

if flowState == "Fully supersonic"
    if s_exit.pressure < atmo.pressure
    flowState = flowState + ", overexpanded";
    else
        flowState = flowState + ", underexpanded";
    end
end

if gas.hasCondensation
    flowState = flowState + ", condensation in the nozzle";
end

displaygasinfo(s_exit, "exit", expansionRatio);

%setsubsonicexitconditions(gas, s_chamber, s_atmo);
specificImpulse_s = specificimpulse(gas, atmo, true);
specificImpulse_ms = specificimpulse(gas, atmo, false);
massFlow = thrust / specificImpulse_ms
A_chamber = massFlow / s_chamber.massFlowFlux; 
heatPower = heatPowerArea * A_chamber;
end

function displaygasinfo(gas, name, areaRatio)
disp(name + ':' + ...
    " Velocity: " + gas.velocity + ...
    " Entropy: " + gas.entropy + ...
    " Total Energy: " + gas.totEnergy + ...
    " Mass Flow Rate: " + gas.massFlowFlux * areaRatio)
end
