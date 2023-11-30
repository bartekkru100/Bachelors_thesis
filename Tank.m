classdef Tank < RocketElement
    properties
        capacity
    end

    methods
        function tank = Tank(prop, diameter, pressure, temperature, capacity)
            tank = tank.defaultconstructor(diameter, diameter)
            tank.inState = State(prop);
            tank.outState = tank.outState
        end
    end
end