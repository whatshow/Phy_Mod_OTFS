classdef OTFS < handle
    properties
        nSubcarNum {mustBeInteger}                  % subcarrier number
        nTimeslotNum {mustBeInteger}                % timeslot number
        freq_spacing {mustBeInteger} = 15           % frequency spacing (kHz), the default is 15kHz
        fc {mustBeInteger} = 3                      % single carrier frequency (GHz), the default is 3GHz
        x_DD                                        % Tx value in the delay Doppler(DD) domain
        x_TF                                        % Tx value in the time-frequency(TF) domain
        s (1, 1) {mustBeVector}                     % Tx value in the time domain (array)
        H                                           % channel in the time domain
        Heff                                        % Effective channel in the DD domain
        r (1, 1) {mustBeVector}                     % Rx value in the time domain (array)
        Y_TF                                        % Rx value in the TF domain
        Y_DD                                        % Rx value in the DD domain
        taps_num                                    % paths number              
        delay_taps (1, 1) {mustBeVector}            % delay index
        dopper_taps (1, 1) {mustBeVector}           % doppler index (it can be integers or fractional numbers)
        chan_coef (1, 1) {mustBeVector}             % path gain
    end
    
    methods
        % constructor
        % @nSubcarNum:      subcarrier number
        % @nTimeslotNum:    timeslot number
        % @freq_spacing:    frequency spacing (kHz), the default is 15kHz
        % @fc:              single carrier frequency (GHz), the default is 3GHz
        function self = OTFS(nSubcarNum, nTimeslotNum, varargin)
            % Inputs Name-Value Pair 
            inPar = inputParser;
            addParameter(inPar,'freq_spacing', self.freq_spacing, @isnumeric);   % register "freq_spacing"
            addParameter(inPar,'fc', self.fc, @isnumeric);                       % register "fc"
            inPar.KeepUnmatched = true;                                          % Allow unmatched cases
            inPar.CaseSensitive = false;                                         % Allow capital or small characters
            
            % take inputs
            % nSubcarNum  
            if isscalar(nSubcarNum) 
                self.nSubcarNum = nSubcarNum;
            else
                error("The number of subcarriers must be a scalar.");
            end
            % nSubcarNum
            if isscalar(nTimeslotNum)
                self.nTimeslotNum = nTimeslotNum;
            else
                error("The number of timeslots must be a scalar.");
            end
            % freq_spacing & fc
            parse(inPar, varargin{:}); 
            self.freq_spacing = inPar.Results.freq_spacing;
            self.fc = inPar.Results.fc;
        end

        % modulate
        % @x_DD: a vector of symbols to send
        function modulate(self, x_DD)
            x_DD = x_DD(:);
            self.x_DD = reshape(x_DD, self.nTimeslotNum, self.nSubcarNum);
            self.x_TF = fft(ifft(self.x_DD).').'/sqrt(M/N); %%%ISFFT
            s_mat = ifft(self.x_TF.')*sqrt(M); % Heisenberg transform
            self.s = s_mat(:);
        end

        % set channel
        % @lmax: the maxmimal delay index
        % @kmax: the maximal Doppler index
        % @p: the path number
        % @is_fractional_doppler: whether we use the fractional Doppler (the default is false)
        function set_channel(self, p, lmax, kmax)
            % Inputs Name-Value Pair 
            inPar = inputParser;
            addParameter(inPar,'is_fractional_doppler', false, @islogical);      % register "is_fractional_doppler"
            inPar.KeepUnmatched = true;                                          % Allow unmatched cases
            inPar.CaseSensitive = false;                                         % Allow capital or small characters

            % load inputs
            if ~isscalar(lmax)
                error("The maximal delay index must be a scalar.");
            elseif lmax >= self.nSubcarNum
                error("The maximal delay index must be less than the subcarrier number.");
            end
            if ~isscalar(kmax)
                error("The maximal Doppler index must be a scalar.");
            elseif kmax > floor(self.nTimeslotNum/2)
                error("The maximal Doppler index must be less than the half of the timeslot number")
            end
            if ~isscalar(p)
                error("The path number must be a scalar.");
            elseif p > lmax*(2*kmax+1)
                error("The path number must be less than lmax*(2*kmax+1) = %d", lmax*(2*kmax+1));
            end
        end
    end

end