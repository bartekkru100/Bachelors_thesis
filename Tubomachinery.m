clc, clear, cleanup

fuel = GRI30;
set(fuel, 'X', 'H2:1')

oxidiser = GRI30;
set(fuel, 'X', 'O2:1')

element = RocketElement(9,9)

$oxidiserTank = Tank(oxidiser, 9, 100e5, 200, 1);
$fuelTank = Tank(fuel, 9, 100e5, 200, 1);