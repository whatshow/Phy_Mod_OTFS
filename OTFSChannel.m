classdef OTFSChannel < handle
    properties
        config;                 % the config of OTFS
    end
    methods
        %{
        constructor
        @config:                the configuration of the entire OTFS sysmtem
        %}
        function self = OTFSChannel(varargin)
            % optional inputs
            % optional inputs - register
            inPar = inputParser;
            addParameter(inPar,"config", self.config, @(x) isscalar(x)&isnumeric(x));
            inPar.KeepUnmatched = true;
            inPar.CaseSensitive = false;
            parse(inPar, varargin{:});
            % optional inputs - load
            self.config = inPar.Results.config;
        end

       
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % obsolete methods
    methods(Static)
        %{
        channel estimation prepare methods
        %}
        function his_mat = obsCSIList2Mat(his, lis, kis, lmax, kmin, kmax)
        end
        function [his, lis, kis] = obsCSIMat2List(his_mat, lmax, kmin, kmax)
        end

        %{
        %}
        function obsCEX(X_DD, lmax, kmin, kmax)
        end
    end
end