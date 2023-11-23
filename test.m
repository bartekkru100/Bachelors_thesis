import State.*
% Test file to test parts of code

clc, clear, cleanup
temperature_in = 293;
pressure_in = 20.64e6;
prop1 = Solution('R1highT.yaml');
set(prop1, 'Y', 'POSF7688:1,O2:2.72', 'P', pressure_in, 'T', temperature_in);
equilibrate(prop1, 'HP');
prop2 = Solution('R1highT.yaml');
equilibrate(prop2, 'HP');
state1 = State(prop1);
state2 = State(prop2);
k(state1);
soundspeed(prop1)