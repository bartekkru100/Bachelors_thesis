classdef State

    % Property names are the same as Cantera's Solution class method names
    % expcept the _mass suffix. State class properties always refer to the
    % mass variant. Velocity property was added for convenience.

    properties (Access = public)
        temperature;
        pressure;
        density;
        meanMolecularWeight;
        enthalpy;
        intEnergy;
        entropy;
        gibbs;
        soundspeed;
        massFractions;
        molecularWeights;
        velocity;
        cp;
        cv;
        k;
        totEnergy;
        Mach;
        R_specific;
        massFlowFlux;

    end
    methods (Access = public)

        % Constructor takes a previously defined Solution or State Object
        % and copies its properties

        function state = State(gas)
            
            % Path for default constructor
            if nargin == 0
                state.temperature = 0;
                state.pressure = 0;
                state.density = 0;
                state.meanMolecularWeight = 0;
                state.enthalpy = 0;
                state.intEnergy = 0;
                state.entropy = 0;
                state.gibbs = 0;
                state.soundspeed = 0;
                state.massFractions = 0;
                state.molecularWeights = 0;
                state.velocity;
                state.cp = 0;
                state.cv = 0;
                state.k = 0;
                state.totEnergy = 0;
                state.Mach = 0;
                state.R_specific = 0;
                state.massFlowFlux = 0;

            % Path for a State class input
            elseif ismember(class(gas), ["Gas", "State"])
                state.temperature = gas.temperature;
                state.pressure = gas.pressure;
                state.density = gas.density;
                state.meanMolecularWeight = gas.meanMolecularWeight;
                state.enthalpy = gas.enthalpy;
                state.intEnergy = gas.intEnergy;
                state.entropy = gas.entropy;
                state.gibbs = gas.gibbs;
                state.soundspeed = gas.soundspeed;
                state.massFractions = gas.massFractions;
                state.molecularWeights = gas.molecularWeights;
                state.velocity = gas.velocity;
                state.cp = gas.cp;
                state.cv = gas.cv;
                state.k = gas.k;
                state.totEnergy = gas.totEnergy;
                state.Mach = gas.Mach;
                state.R_specific = gas.R_specific;
                state.massFlowFlux = gas.massFlowFlux;
            end
        end
    end
end