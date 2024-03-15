classdef MullersMethod < RootFindingMethod

    % Muller's method

    properties (Access = private)
        min_max;
    end

    methods (Access = public)

        function this = MullersMethod(location, tolerance, iterationLimit, X, Y, min_max)
            this@RootFindingMethod(location, tolerance, iterationLimit, X, Y);
            this.min_max = Enum('min', 'max');
            this.min_max.current = this.min_max.fromstr(min_max);
        end

    end

    methods (Access = protected)

        function X_0 = findnewX_subclass(this)
            p = this.X(1) - this.X(2);
            q = p / (this.X(2) - this.X(3));
            a = q * (this.Y(1) - (1 + q) * this.Y(2)) + q ^ 2 * this.Y(3);
            b = (2 * q + 1) * this.Y(1) - (1 + q) ^ 2 * this.Y(2) + q ^ 2 * this.Y(3);
            c = (1 + q) * this.Y(1);
            sqrtDelta = sqrt(b ^ 2 - 4 * a * c);
            denom(1) = (b + sqrtDelta);
            denom(2) = (b - sqrtDelta);
            if this.min_max.current == this.min_max.value.max
                X_0 = max(this.X(1) - p * (2 * c) ./ denom);
            elseif this.min_max.current == this.min_max.value.min
                X_0 = min(this.X(1) - p * (2 * c) ./ denom);
            else
                error("Specify 'min' or 'max'");
            end
        end

        function updateXY_subclass(this, X_0, Y_0)
            this.X(3) = this.X(2);
            this.Y(3) = this.Y(2);
            this.X(2) = this.X(1);
            this.Y(2) = this.Y(1);
            this.X(1) = X_0;
            this.Y(1) = Y_0;
        end

    end
end