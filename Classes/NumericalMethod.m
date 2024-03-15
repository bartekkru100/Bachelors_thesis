classdef NumericalMethod < handle

    % This is a virtual class used as a template for root finding methods

    properties (Access = public)
        location; % An identifier to make it easier to see where convergence failed, helps with debugging.
        tolerance; % Required tolerance.
        iterationLimit; % Self explanatory.
        nPoints; % Number of initial points required to get a root approximation.
        convergenceCounter; % Keeps tack of whether or not the algorithm is converging.
        shouldBreak; % Boolean decides whether or not the algorithm should terminate when checkconvergence() is called.
        warningsEnabled; % Flag indicates whether or not warnings should be displayed on non convergence.
        hasFailed; % Flag indicates whether or not convergence failed.

        % Current root approximation
        Y_0;
        X_0;
        error_0;

        % Approximation from the previous iteration
        Y_0_old;
        X_0_old;
        error_0_old;

        % Approximation bounds, any points calculated beyond those will
        % need to be recalculated using a different method
        X_min
        X_max;

        % Points used for calculating the root
        X;
        Y;
        
        nIterations; % Number of iterations
    end

    methods (Access = public)

%==========================================================================

        % Class constructor

        function this = NumericalMethod(location, tolerance, iterationLimit, X, Y)
            import NumericalMethod.*
            this.location = location;
            this.tolerance = tolerance;
            this.iterationLimit = iterationLimit;
            this.error_0_old = tolerance * 100;
            this.error_0 = tolerance * 10;
            this.nIterations = 0;
            this.convergenceCounter = 1;
            this.shouldBreak = false;
            this.hasFailed = false;
            this.warningsEnabled = true;
            if nargin == 5
                this.X = X;
                this.Y = Y;
            end
        end

%==========================================================================

        % Checks for convergence and eventually breaks the loop

        function shouldBreak = checkconvergence(this)
            this.nIterations = this.nIterations + 1;
            this.errorCalculation;
            if this.shouldBreak
            elseif abs(this.error_0) < this.tolerance
                this.shouldBreak = true;
            elseif this.nIterations == 1
            else
                convergence = abs(2 * this.error_0 / this.error_0_old);
                if this.nIterations > this.iterationLimit
                    this.hasFailed = true;
                elseif convergence >= 1
                    if this.convergenceCounter < 1
                        this.convergenceCounter = 1;
                    end
                    this.convergenceCounter = this.convergenceCounter * convergence;
                    if this.convergenceCounter > 1000
                        this.hasFailed = true;
                    end
                elseif convergence < 1
                    this.convergenceCounter = this.convergenceCounter * sqrt(convergence);
                end
            end
            if this.hasFailed
                if this.warningsEnabled
                    printconvergencewarning(this, this.error_0);
                end
                this.shouldBreak = true;
            end
            shouldBreak = this.shouldBreak;
        end

%==========================================================================
        
        % Checks whether convergence happened earlier

        function hasFailed = hasfailed(this)
            hasFailed = this.hasFailed;
        end

%==========================================================================

        % Updates X and Y using the new calculated values

        function updateXY(this, X_0, Y_0)
            this.X_0 = X_0;
            this.Y_0 = Y_0;
            updateXY_subclass(this, X_0, Y_0); % Call to a virtual method defined in a subclass
        end

%==========================================================================

        % Gives or takes a new approximation of the root. 

        function X_0 = findnewX(this, varargin)
            this.Y_0_old = this.Y_0;
            this.X_0_old = this.X_0;
             % Calls to a virtual method defined in a subclass
            if nargin == 1
                this.X_0 = real(findnewX_subclass(this));
            else
                this.X_0 = real(findnewX_subclass(this, varargin));
            end

%--------------------------------------------------------------------------

            % This potentially applies limits to the newly calculated root

            if this.X_0 < this.X_min
                if ~isempty(this.X)
                    this.X_0 = (min(this.X) + this.X_min) / 2;
                elseif ~isempty(this.X_0_old)
                    this.X_0 = (this.X_0_old + this.X_min) / 2;
                else
                    this.X_0 = (this.X_max + 3 * this.X_min) / 4;
                end
            elseif this.X_0 > this.X_max
                if ~isempty(this.X)
                    this.X_0 = (max(this.X) + this.X_max) / 2;
                elseif ~isempty(this.X_0_old)
                    this.X_0 = (this.X_0_old + this.X_max) / 2;
                else
                    this.X_0 = (3 * this.X_max + this.X_min) / 4;
                end
            elseif isnan(this.X_0)
                %this.X_0 = this.X_0_old;
                this.X_0 = mean(this.X);
                %this.hasFailed = true;
                %this.shouldBreak = true;
            end
            X_0 = this.X_0;
        end

%==========================================================================

        % Some getters and setters

        function setX_min(this, X_min)
                this.X_min = X_min;
        end

        function setX_max(this, X_max)
                this.X_max = X_max;
        end

        function setX_min_max(this, X_min, X_max)
                this.X_min = X_min;
                this.X_max = X_max;
        end

        function XY = get(this)
            XY.X = this.X;
            XY.Y = this.Y;
            XY.X_0 = this.X_0;
            XY.Y_0 = this.Y_0;
            XY.X_0_old = this.X_0_old;
            XY.Y_0_old = this.Y_0_old;
        end

        function setnewX(this, X_0)
            this.X_0 = X_0;
        end

        function iteration = iteration(this)
            iteration = this.nIterations + 1;
        end

        function disablewarnings(this)
            this.warningsEnabled = false;
        end

        function enablewarnings(this)
            this.warningsEnabled = true;
        end

%==========================================================================

        % WIP, to be remade. This will faciliate switching between numerical methods on
        % the fly

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

%==========================================================================

    % Virtual methods called inside the main class methods

    methods (Access = protected, Abstract)

        errorCalculation(this)

        findnewX_subclass(this, varargin)

        updateXY_subclass(this, X_0, Y_0)

    end

    methods(Access = protected)

%==========================================================================

        % Prints a warning if convergence fails

        function printconvergencewarning(this, error)
            warning("convergence failed at " + this.location + " at iteration " + this.nIterations + " with error = " + error);
        end

    end
end