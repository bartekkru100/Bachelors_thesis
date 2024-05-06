function solverOutputs = solvererroroutput
global solverError

solverOutputs.effectiveExhaustVelocity = [];
solverOutputs.specificImpulse = [];
solverOutputs.thrustCoefficient = [];
solverOutputs.characteristicVelocity = [];
solverOutputs.A_chamber = [];
solverOutputs.A_throat = [];
solverOutputs.A_exit = [];
solverOutputs.D_chamber = [];
solverOutputs.D_throat = [];
solverOutputs.D_exit = [];
solverOutputs.massFlow = [];
solverOutputs.heatPower = [];
solverOutputs.shockPosition = [];
solverOutputs.s_shock_1 = [];
solverOutputs.s_shock_2 = [];
solverOutputs.s_injection = [];
solverOutputs.s_chamber = [];
solverOutputs.s_throat = [];
solverOutputs.s_exit = [];
solverOutputs.optimalRatio = [];
solverOutputs.pressure_ideal = [];
solverOutputs.pressure_shockAtExit = [];
solverOutputs.pressure_shockAtThroat = [];
solverOutputs.pressure_separation = [];
solverOutputs.hasCondensation = [];
solverOutputs.expansionRatio_shockAtThroat = [];
solverOutputs.expansionRatio_optimal = [];
solverOutputs.expansionRatio_shockAtExit_1 = [];
solverOutputs.expansionRatio_shockAtExit_2 = [];
solverOutputs.expansionRatio_condensation = [];
solverOutputs.expansionRatio_separation = [];
solverOutputs.flowState = solverError;
end