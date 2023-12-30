clc, clear

Table = readtable("SteamTable.txt");
temperature = table2array(Table(:,1));
pressure = table2array(Table(:,2));
temps = linspace(273.15, 500, 20);
press = interp1(temperature, pressure, temps, "spline");
plot(temps, press)