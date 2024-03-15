classdef LocalMinMax < NumericalMethod

    % Muller's method

    properties (Access = private)
        min_max;
        up_down;
    end

    methods (Access = public)

        function this = LocalMinMax(location, tolerance, iterationLimit, X, Y, min_max)
            this@NumericalMethod(location, tolerance, iterationLimit, X, Y);
            this.min_max = Enum('min', 'max');
            this.min_max.current = this.min_max.fromstr(min_max);
            this.up_down = Enum('up', 'down');
        end

    end

    methods (Access = protected)

        function errorCalculation(this)
            this.error_0_old = this.error_0;
            if this.min_max.current == this.min_max.fromstr('min')
                this.error_0 = (min(this.Y) - max(this.Y)) / min(this.Y);
            else
                this.error_0 = (max(this.Y) - min(this.Y)) / max(this.Y);
            end
        end

        function X_0 = findnewX_subclass(this)
            p = this.X(2) - this.X(1);
            q = this.X(2) ^ 2 - this.X(1) ^ 2;
            r = (this.X(3) - this.X(1)) / p;
            a = (this.Y(3) - this.Y(1) - (this.Y(2) - this.Y(1)) * r) / ...
                (this.X(3) ^ 2 - this.X(1) ^ 2 - q * r);
            if a > 0
                this.up_down.current = this.up_down.fromstr('up');
            else
                this.up_down.current = this.up_down.fromstr('down');
            end

            if this.min_max.current == this.up_down.current
                X_0 = - (this.Y(2) - this.Y(1) - a * q) / (2 * a * p);
            else
                if this.min_max.current == this.min_max.fromstr('min')
                    [Y_sorted, I] = sort(this.Y, 'ascend');
                elseif this.min_max.current == this.min_max.fromstr('max')
                    [Y_sorted, I] = sort(this.Y, 'descend');
                end
                X_sorted = this.X(I);
                X_0 = 3 * X_sorted(1) - 2 * X_sorted(2);
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