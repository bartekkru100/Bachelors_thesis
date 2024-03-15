classdef ErrorCorrectedMethod < RootFindingMethod

    % Error corrected method

    properties (Access = protected)
        Y_abs;
        Y_abs_old;
    end

    methods (Access = protected)

        function X_0 = findnewX_subclass(this, X_0)
            X_0 = cell2mat(X_0);
        end

        function updateXY_subclass(this, X_0, Y_0)
            this.X = X_0;
            this.Y = Y_0;
            if nargout == 2
                X = this.X;
                Y = this.Y;
            end
        end

    end
end