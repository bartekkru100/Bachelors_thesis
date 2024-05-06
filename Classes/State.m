classdef State

    % Property names are the same as Cantera's Solution class method names
    % expcept the _mass suffix. State class properties always refer to the
    % mass variant. Velocity property was added for convenience.

    properties (Access = public)
        src;
        temperature;
        pressure;
        density;
        enthalpy;
        intEnergy;
        entropy;
        gibbs;
        totEnergy;
        velocity;
        Mach;
        soundspeed;
        massFlowFlux;
        areaRatio
        cp;
        cv;
        k;
        R_specific;
        meanMolecularWeight;
        hasCondensation;
        massFractions;
    end
    methods (Access = public)

        % Constructor takes a previously defined Solution or State Object
        % and copies its properties

        function state = State(gas)
            
            % Path for default constructor
            if nargin == 0
                state.src = "GRI30";
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
                state.velocity = 0;
                state.cp = 0;
                state.cv = 0;
                state.k = 0;
                state.totEnergy = 0;
                state.Mach = 0;
                state.R_specific = 0;
                state.massFlowFlux = 0;
                state.areaRatio = 0;
                state.hasCondensation = false;

            % Path for a State class input
            elseif ismember(class(gas), ["Gas", "State"])
                state.src = gas.src;
                state.temperature = gas.temperature;
                state.pressure = gas.pressure;
                state.density = gas.density;
                state.meanMolecularWeight = gas.meanMolecularWeight;
                state.enthalpy = gas.enthalpy;
                state.intEnergy = gas.intEnergy;
                state.entropy = gas.entropy;
                state.gibbs = gas.gibbs;
                state.totEnergy = gas.totEnergy;
                state.velocity = gas.velocity;
                state.Mach = gas.Mach;
                state.soundspeed = gas.soundspeed;
                state.massFlowFlux = gas.massFlowFlux;
                state.areaRatio = gas.areaRatio;
                state.cp = gas.cp;
                state.cv = gas.cv;
                state.k = gas.k;
                state.R_specific = gas.R_specific;
                state.massFractions = gas.massFractions;
                state.hasCondensation = gas.hasCondensation;
            end
        end
    end
end