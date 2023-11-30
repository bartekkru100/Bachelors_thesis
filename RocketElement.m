classdef RocketElement
    properties
        name;
        inDiameter;
        outDiameter;
        inState;
        outState;
        massFlow;
    end

    methods
        function element = RocketElement(inDiameter,outDiameter)
            element = element.defaultconstructor(element, inDiameter, outDiameter)
        end

        function element = defaultconstructor(element, inDiameter, outDiameter)
            element.inDiameter = inDiameter;
            element.outDiameter = outDiameter;
        end
    end
end