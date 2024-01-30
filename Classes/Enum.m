classdef Enum

    % This class is used as a replacement for Matlab's own enumerations 

    properties (Access = protected)
        values;
        names;
    end

    properties (Access = public)
        current;
    end
    
    methods (Access = public)
        
        function this = Enum(varargin)
            this.names = varargin;
            for i = 1 : size(varargin, 2)
                this.values.(cell2mat(varargin(i))) = i;
            end
        end

        function value = value(this)
            value = this.values;
        end

        function value = fromstr(this, str)
            str = convertStringsToChars(str);
            value = this.values.(str);
        end

        function str = tostr(this, value)
            str = cell2mat(this.names(value));
        end

    end

end