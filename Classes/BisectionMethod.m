classdef BisectionMethod < RootFindingMethod

    % Bisection method

    methods (Access = protected)

        function X_0 = findnewX_subclass(this)
            X_0 = (this.X(1) + this.X(2)) / 2;
        end

        function updateXY_subclass(this, X_0, Y_0)

            % New points are chosen based on the sign of error  

            if sign(Y_0) == sign(this.Y(1))
                this.X(1) = X_0;
                this.Y(1) = Y_0;
            else
                this.X(2) = X_0;
                this.Y(2) = Y_0;
            end
            if nargout == 2
                X = this.X;
                Y = this.Y;
            end
        end

    end
end