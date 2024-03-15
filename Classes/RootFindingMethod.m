classdef RootFindingMethod < NumericalMethod

    % Bisection method

    methods (Access = protected)

        function errorCalculation(this)
            this.error_0_old = this.error_0;
            this.error_0 = this.Y_0;
        end

    end

%==========================================================================

    % Virtual methods called inside the main class methods

    methods (Access = protected, Abstract)

        findnewX_subclass(this, varargin)

        updateXY_subclass(this, X_0, Y_0)

    end
end