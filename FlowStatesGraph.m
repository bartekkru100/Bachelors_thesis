addpathall(pwd, ".git", true);
addpath('D:\Uni\Inzynierka\Program')
clc, clear%, cleanup
import Gas.*

atmo = Gas('air');
gas = Gas('GRI30');
i = 0;
for expansionRatio = logspace(0,log10(3),100)
    i = i + 1;
    setstate(gas, 'Y', 'H2:1', 'T', 50, 'P', 20e5, 'velocity', 0);
    propellantArray(1) = State(gas);
    setstate(gas, 'Y', 'O2:1', 'T', 50, 'P', 20e5, 'velocity', 0);
    propellantArray(2) = State(gas);
    setstate(atmo, 'P', oneatm, 'T', atmo.temperature);

    solverInputs.atmo = atmo;
    solverInputs.gas = gas;
    solverInputs.propellantArray = propellantArray;
    solverInputs.ratios = [1, 6];
    solverInputs.contractionRatio = 5;
    solverInputs.expansionRatio = expansionRatio;
    solverInputs.thrust = 1.0331e6;
    solverInputs.heatMass = 10000000;
    solverInputs.pressure_injection = 200e5;
    solverInputs.temperature_injection = 50;
    solverInputs.temperature_chamber_max = 4000;
    solverInputs.temperature_throat_max = 4000;
    solverInputs.separationTolerance = 0;
    solverInputs.mass_mole = "mass";
    solverInputs.searchedIndicator = "specificImpulse";
    solverOutputs = engineanalysis(solverInputs);
    % solverOutputs = optimalnozzle(solverInputs);
    % solverOutputs = optimalmixture(solverInputs);
    x(i) = expansionRatio;
    y(i, 1) = solverOutputs.pressure_shockAtThroat / solverInputs.pressure_injection;
    y(i, 2) = solverOutputs.pressure_shockAtExit / solverInputs.pressure_injection;
    y(i, 3) = solverOutputs.pressure_separation / solverInputs.pressure_injection;
    y(i, 4) = solverOutputs.pressure_ideal / solverInputs.pressure_injection;
end
y = y;

fig = figure('Position', [0, 0, 600, 800]);
tiledlayout(1,1, 'Padding', 'compact');

nexttile
plot(x, y(:, 1), "Color", "black");
hold on
plot(x, y(:, 2), "Color", "black");
hold on
plot(x, y(:, 3), "Color", "black", "LineStyle","--");
hold on
plot(x, y(:, 4), "Color", "black");
ylim padded

yline(gas.sonic.pressure / solverInputs.pressure_injection, ...
    "--", "$P_{throat}$", "Color", "black", "Interpreter", "latex", ...
    "LabelVerticalAlignment", "middle", "LabelHorizontalAlignment", "left");

yline(gas.stagnation.pressure / solverInputs.pressure_injection, ...
    "-", "$P_{0}$", "LineWidth", 1, "Color", "black", "Interpreter", "latex", ...
    "LabelVerticalAlignment", "middle", "LabelHorizontalAlignment", "left");

text(2, 1.03, ...
    "Flow impossible", ...
    "VerticalAlignment", "middle", "HorizontalAlignment", "center", ...
    "Rotation", 0);

text(x(20), mean([y(20,1), 1]), ...
    "Subsonic", ...
    "VerticalAlignment", "middle", "HorizontalAlignment", "center", ...
    "Rotation", 0);

text(x(75), mean([y(75,2), y(75,1)]),  ...
    "Normal shock in the nozzle", ...
    "VerticalAlignment", "middle", "HorizontalAlignment", "center", ...
    "Rotation", 0);

text(x(60), mean([y(60,3), y(60,2)]), ...
    "Flow separation", ...
    "VerticalAlignment", "middle", "HorizontalAlignment", "center", ...
    "Rotation", 0);

text(x(70), mean([y(70,4), y(70,3)]), ...
    "Overexpanded", ...
    "VerticalAlignment", "middle", "HorizontalAlignment", "center", ...
    "Rotation", 0);

text(x(35), mean([0, y(35,4)]), ...
    "Underexpanded", ...
    "VerticalAlignment", "middle", "HorizontalAlignment", "center", ...
    "Rotation", 0);

text(x(50), y(50,1), ...
    "Shock at throat", ...
    "VerticalAlignment", "middle", "HorizontalAlignment", "center", ...
    "BackgroundColor", "white", "Rotation", 16);

text(x(30), y(30,2), ...
    "Shock at exit", ...
    "VerticalAlignment", "middle", "HorizontalAlignment", "center", ...
    "BackgroundColor", "white", "Rotation", -26);

text(x(40), y(40,3), ...
    "Separation threshold", ...
    "VerticalAlignment", "middle", "HorizontalAlignment", "center", ...
    "BackgroundColor", "white", "Rotation", -32);

text(x(40), y(40,4), ...
    "Ideal expansion", ...
    "VerticalAlignment", "middle", "HorizontalAlignment", "center", ...
    "BackgroundColor", "white", "Rotation", -29);

xlabel('$A_e/{A^*}$', 'Interpreter', 'latex');
ylabel("$P_{atm}/P_0$", 'Interpreter', 'latex');

title("Flow regimes", 'Interpreter', 'latex')

fontname("Cambria Math")
fontsize(11, 'points')
