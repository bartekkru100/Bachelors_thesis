classdef State
    properties (Access = public)
        temperature;
        pressure;
        density;
        meanMolecularWeight;
        enthalpy;
        intEnergy;
        entropy;
        gibbs;
        cp;
        cv;
        soundspeed;
        massFractions;
        molecularWeights;
        velocity;

    end
    methods (Access = public)
        function state = State(gas, velocity)
            if class(gas) == "Solution"
                state.temperature = temperature(gas);
                state.pressure = pressure(gas);
                state.density = density(gas);
                state.meanMolecularWeight = meanMolecularWeight(gas);
                state.enthalpy = enthalpy_mass(gas);
                state.intEnergy = intEnergy_mass(gas);
                state.entropy = entropy_mass(gas);
                state.gibbs = gibbs_mass(gas);
                state.cp = cp_mass(gas);
                state.cv = cv_mass(gas);
                state.soundspeed = soundspeed(gas);
                state.massFractions = massFractions(gas);
                state.molecularWeights = molecularWeights(gas);
                if nargin == 2
                    state.velocity = velocity;
                end
            elseif class(gas) == "State"
                

            else

            end
        end
    end
    methods (Static)
        function setState(gas, state)
            set(gas, 'Y', state.massFractions, 'T', state.temperature,'P', state.pressure);
        end
    end
end