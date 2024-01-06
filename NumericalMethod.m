classdef NumericalMethod < handle

    properties (Access = protected)
        location;
        tolerance;
        iterationLimit;
        nIterations;
        nPoints;
        convergenceCounter;

        Y_0_old;
        X_0_old;

        Y_0;
        X_0;

        X;
        Y;
    end

    methods (Access = public)

        function this = NumericalMethod(location, tolerance, iterationLimit, X, Y)
            import NumericalMethod.*

            this.location = location;
            this.tolerance = tolerance;
            this.iterationLimit = iterationLimit;
            this.Y_0_old = tolerance * 10;
            this.nIterations = 0;
            this.convergenceCounter = 0;
            if nargin == 5
                this.X = X;
                this.Y = Y;
            end
        end

        function shouldBreak = checkconvergence(this)
            this.nIterations = this.nIterations + 1;
            shouldBreak = false;
            if abs(this.Y_0) < this.tolerance
                shouldBreak = true;
            elseif this.nIterations > this.iterationLimit
                printconvergencewarning(this, this.Y_0);
                shouldBreak = true;
            elseif this.nIterations == 1

            elseif abs(this.Y_0 / this.Y_0_old) >= 1
                if this.convergenceCounter < 0
                    this.convergenceCounter = 1;
                else
                    this.convergenceCounter = this.convergenceCounter + 2;
                end
                if this.convergenceCounter > 7
                    printconvergencewarning(this, this.Y_0);
                    shouldBreak = true;
                end
            else
                this.convergenceCounter = this.convergenceCounter - 1;
            end
        end

        function updateXY(this, X_0, Y_0)
            this.X_0 = X_0;
            this.Y_0 = Y_0;
            updateXY_subclass(this, X_0, Y_0);
        end

        function X_0 = findnewX(this, varargin)
            this.Y_0_old = this.Y_0;
            this.X_0_old = this.X_0;
            if nargin == 1
                this.X_0 = findnewX_subclass(this);
            else
                this.X_0 = findnewX_subclass(this, varargin);
            end
            X_0 = this.X_0;
        end

        function hasConverged = hasconverged(this)
            hasConverged = this.Y_0 < this.tolerance;
        end

        function XY = get(this)
            XY.X = this.X;
            XY.Y = this.Y;
            XY.X_0 = this.X_0;
            XY.Y_0 = this.Y_0;
            XY.X_0_old = this.X_0_old;
            XY.Y_0_old = this.Y_0_old;
        end

        function iteration = iteration(this)
            iteration = this.nIterations + 1;
        end

        function reset(this, method, X, Y)
            method = toMethod(method);
            this.nIterations = 0;
            this.convergenceCounter = 0;
            switch numberofpoints(method)
                case 1
                    if nargin == 1
                        this.X = this.X(1);
                        this.Y = this.Y(1);
                    elseif nargin == 3
                        this.X = X(1);
                        this.Y = Y(1);
                    else
                        error("Incorrect number of points for: " + methodToString(method))
                    end
                case 2
                    if nargin == 1 & size(this.X, 2) == 2
                    elseif nargin == 1 & size(this.X, 2) < 2
                        [this.X, order] = sort(this.X);
                        this.Y = this.Y(order);
                        this.X(3 : size(this.X, 2)) = [];
                        this.Y(3 : size(this.Y, 2)) = [];
                    elseif nargin == 3 & size(X, 2) == 2
                        this.X = X;
                        this.Y = Y;
                    elseif nargin == 3 & size(X, 2) == 1 & size(this.X, 2) == 1
                        this.X(2) = X;
                        this.Y(2) = Y;
                    else
                        error("Incorrect number of points for: " + methodToString(method))
                    end
                case 3
                    if nargin == 1 & size(this.X, 2) == 3
                    elseif nargin == 1 & size(this.X, 2) < 3
                        [this.X, order] = sort(this.X);
                        this.Y = this.Y(order);
                        this.X(4 : size(this.X, 2)) = [];
                        this.Y(4 : size(this.Y, 2)) = [];
                    elseif nargin == 1 & size(this.X, 2) > 3
                        this.X(3) = this.X(2);
                        this.Y(3) = this.Y(2);
                        this.X(2) = (this.X(1) + this.X(3)) / 2;
                        this.Y(2) = (this.X(1) + this.X(3)) / 2;
                    elseif nargin == 3 & size(X, 2) == 3
                        this.X = X;
                        this.Y = Y;
                    elseif nargin == 3 & size(X, 2) == 2 & size(this.X, 2) == 1
                        this.X(3) = X;
                        this.Y(3) = Y;
                    elseif nargin == 3 & size(X, 2) == 1 & size(this.X, 2) == 2
                        this.X(2 : 3) = X;
                        this.Y(2 : 3) = Y;
                    else
                        error("Incorrect number of points for: " + methodToString(method))
                    end
            end
            this.method = method;
        end

    end

    methods (Access = protected, Abstract)

        findnewX_subclass(this, varargin)

        updateXY_subclass(this, X_0, Y_0)

    end

    methods(Access = protected)

        function printconvergencewarning(this, error)
            warning("convergence failed at " + this.location + " at iteration " + this.nIterations + " with error = " + error);
        end

    end

    methods (Static, Access = private)

        function method = numericalmethods()
            method = struct('errorCorrected', 1, 'bisection', 2, 'Mullers', 3);
        end

        function method = toMethod(str)
            switch lower(convertStringsToChars(lower(str)))
                case 'errorcorrected'
                    method = 1;
                case 'bisection'
                    method = 2;
                case 'mullers'
                    method = 3;
                otherwise
                    error("Wrong input: " + str);
            end
        end

        function str = methodToString(method)
            switch method
                case 1
                    str = 'errorCorrected';
                case 2
                    str = 'bisection';
                case 3
                    str = 'mullers';
                otherwise
                    error("Wrong numerical method number: " + method);
            end
        end

        function nX = numberofpoints(method)
            switch method
                case numericalmethods().errorCorrected
                    nX = 1;
                case numericalmethods().bisection
                    nX = 2;
                case numericalmethods().Mullers
                    nX = 3;
                otherwise
                    error("Wrong numerical method name: " + str);
            end
        end

    end
end