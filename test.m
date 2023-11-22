import State.*
% Test file to test parts of code

clc, clear, cleanup
temperature_in = 293;
pressure_in = 20.64e6;
prop = Solution('R1highT.yaml');
set(prop, 'Y', 'POSF7688:1,O2:2.72', 'P', pressure_in, 'T', temperature_in);
equilibrate(prop, 'HP');

prop_cc = State(prop, 1000)
prop_1 = Solution('R1highT.yaml');
setState(prop_1, prop_cc)