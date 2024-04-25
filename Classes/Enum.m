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
            this.names = string.empty();
            for i = 1 : size(varargin, 2)
            name = lower(varargin{i});
            name = strrep(name, ' ', '_');
            name = strrep(name, ',', '');
            this.names(i) = varargin{i};
            this.values.(name) = i;
            end
        end

        function value = value(this)
            value = this.values;
        end

        function value = fromstr(this, str)
            str = convertCharsToStrings(str);
            value = this.values.(str);
        end

        function str = tostr(this, value)
            str = this.names(value);
        end

    end

end