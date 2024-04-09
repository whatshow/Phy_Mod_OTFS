classdef OTFS < handle
    % constants
    properties(Constant)
        % Pulse Types
        PULSE_IDEAL = 10;   % ideal pulses (if we use ideal pulses, `ISI_CANCEL_CP` is forced to chosen)
        PULSE_RECTA = 20;   % rectangular pulses
        PULSE_TYPES = [OTFS.PULSE_IDEAL, OTFS.PULSE_RECTA];
        % isi cancel methods
        ISI_CANCEL_GUARD = 0;                       % use zero guards to avoid ISI
        ISI_CANCEL_CP = 10;                         % one cp for entire OTFS frame
        ISI_CANCEL_CP_MULTI = 11;                   % one cp for each OTFS subframe
        ISI_CANCEL_ZP = 20;                         % zero padding
        ISI_CANCEL_TYPES = [OTFS.ISI_CANCEL_GUARD, OTFS.ISI_CANCEL_CP, OTFS.ISI_CANCEL_CP_MULTI, OTFS.ISI_CANCEL_ZP];
        % Detect
        DETECT_MP_BASE = 10;                        % base OTFS MP detector proposed by P. Raviteja in 2018
        DETECT_TYPES = [OTFS.DETECT_MP_BASE];
        % CSI
        DETECT_CSI_PERFECT = 10;                    % perfect CSI
        DETECT_CSI_CE = 20;                         % CSI from channel estimation (its type is based on pilot type)
        DETECT_CSI_IN = 30;                         % CSI from the input
        DETECT_CSI_TYPES = [OTFS.DETECT_CSI_PERFECT, OTFS.DETECT_CSI_CE, OTFS.DETECT_CSI_IN];
        % CE Area Tag
        CE_AREA_TAG_XDD = 10;
        CE_AREA_TAG_YDD = 11;
        CE_AREA_TAGS = [OTFS.CE_AREA_TAG_XDD, OTFS.CE_AREA_TAG_YDD];
    end
    properties
        nSubcarNum {mustBeInteger}                  % subcarrier number
        nTimeslotNum {mustBeInteger}                % timeslot number
        X_DD                                        % Tx value in the delay Doppler(DD) domain
        X_TF                                        % Tx value in the time-frequency(TF) domain
        s                                           % Tx value in the time domain (array)
        H                                           % channel in the time domain
        r                                           % Rx value in the time domain (array)
        Y_TF                                        % Rx value in the TF domain
        Y_DD                                        % Rx value in the DD domain
        y_DD                                        % Rx value in the DD domain (vectorized)
        taps_num = 0                                % paths number              
        delay_taps                                  % delay index, a row vector
        doppler_taps                                % doppler index (integers or fractional numbers), a row vector
        chan_coef                                   % path gain, a row vector
        % pulse & ISI
        pulse_type = OTFS.PULSE_RECTA;
        isi_cancel_type = OTFS.ISI_CANCEL_CP;
        % detection
        constellation                               % constellation
        constellation_len = 0;
    end
    
    % General OTFS Methods
    methods   
        %{
        constructor
        @nSubcarNum:                subcarrier number
        @nTimeslotNum:              timeslot number
        @pulse_type:                pulse shaping
        @isi_cancel_type:
        %}
        function self = OTFS(nSubcarNum, nTimeslotNum, varargin)
            % Inputs Name-Value Pair 
            inPar = inputParser;
            addParameter(inPar, "pulse_type", self.pulse_type, @(x) isscalar(x)&&isnumeric(x)&&ismember(x, self.PULSE_TYPES));
            addParameter(inPar, "isi_cancel_type", self.isi_cancel_type, @(x) isscalar(x)&&isnumeric(x)&&ismember(x, self.ISI_CANCEL_TYPES));
            inPar.KeepUnmatched = true;                                             % Allow unmatched cases
            inPar.CaseSensitive = false;                                            % Allow capital or small characters
            parse(inPar, varargin{:});

            % take inputs
            % nSubcarNum  
            if ~isscalar(nSubcarNum) || nSubcarNum ~= floor(nSubcarNum)
                error("The number of subcarriers must be an integer scalar.");
            else
                self.nSubcarNum = nSubcarNum;
            end
            % nSubcarNum
            if ~isscalar(nTimeslotNum) || nTimeslotNum ~= floor(nTimeslotNum)
                error("The number of timeslots must be an integer scalar.");
            else
                self.nTimeslotNum = nTimeslotNum;
            end
            % pulse shaping
            self.pulse_type = inPar.Results.pulse_type;
            % isi cancel type
            self.isi_cancel_type = inPar.Results.isi_cancel_type;
            if self.pulse_type == self.PULSE_IDEAL
                self.isi_cancel_type = self.ISI_CANCEL_CP;
            end
        end
        
        %{
        modulate
        @rg: an OTFS resource grid
        %}
        function modulate(self, rg)
            % ISFFT 
            X_FT = fft(ifft(rg.getXDD()).').'/sqrt(self.nSubcarNum/self.nTimeslotNum);
            % X_TF is [nSubcarNum, nTimeslotNum]
            self.X_TF = X_FT.'; 
            % Heisenberg transform
            s_mat = ifft(self.X_TF)*sqrt(self.nSubcarNum);
            % vectorize
            self.s = s_mat(:);
        end

        %{
        demodulate
        %}
        function yDD = demodulate(self)
            r_mat = reshape(self.r, self.nSubcarNum, self.nTimeslotNum);
            % Wigner transform (Y_TF in [nSubcarNum, nTimeslotNum])
            self.Y_TF = fft(r_mat)/sqrt(self.nSubcarNum); 
            Y_FT = self.Y_TF.';
            % SFFT (Y_DD in [Doppler, delay] or [nTimeslotNum ,nSubcarNum])
            self.Y_DD = ifft(fft(Y_FT).').'/sqrt(self.nTimeslotNum/self.nSubcarNum); 
            % vectorize
            if ~self.isInsertPilotsAndGuards()
                % no invalid area
                yDD = self.Y_DD.';
                yDD = yDD(:);
            else
                % return the invalid
                yDD = zeros(self.nTimeslotNum*self.nSubcarNum - (self.ce_ydd_doppl_end-self.ce_ydd_doppl_beg+1)*(self.ce_ydd_delay_end-self.ce_ydd_delay_beg+1), 1);
                symbols_id = 1;
                for doppl_id = 1:self.nTimeslotNum
                    for delay_id = 1:self.nSubcarNum
                        if doppl_id<self.ce_ydd_doppl_beg || doppl_id>self.ce_ydd_doppl_end || delay_id<self.ce_ydd_delay_beg || delay_id>self.ce_ydd_delay_end
                            yDD(symbols_id) = self.Y_DD(doppl_id, delay_id);
                            symbols_id = symbols_id + 1;
                        end
                    end
                end
            end
            % store
            self.y_DD = yDD;
        end

        % set channel
        % <set random paths>
        % @p: the path number
        % @lmax: the maxmimal delay index
        % @kmax: the maximal Doppler index (can be fractional)
        % <set known paths>
        % @delays:      the delays
        % @dopplers:    the doppler shifts
        % @gains:       the path gains
        function setChannel(self, varargin)
            % Inputs Name-Value Pair 
            inPar = inputParser;
            addParameter(inPar,'p', 0, @(x) isscalar(x)&&round(x)==x);
            addParameter(inPar,'lmax', 0, @(x) isscalar(x)&&round(x)==x);
            addParameter(inPar,'kmax', 0, @isscalar);
            addParameter(inPar,'delays', [], @isnumeric);
            addParameter(inPar,'dopplers', [], @isnumeric);
            addParameter(inPar,'gains', [], @isnumeric);
            inPar.KeepUnmatched = true;                                          % Allow unmatched cases
            inPar.CaseSensitive = false;                                         % Allow capital or small characters
            % load inputs
            parse(inPar, varargin{:});
            p = inPar.Results.p;
            lmax = inPar.Results.lmax;
            kmax = inPar.Results.kmax;
            kmax_frac = 0;
            is_fractional_doppler = false;
            if kmax ~= floor(kmax)
                is_fractional_doppler = true;
                kmax_frac = kmax;
                kmax = floor(kmax);
            end
            delays = inPar.Results.delays;
            dopplers = inPar.Results.dopplers;
            gains = inPar.Results.gains;
            
            % check which input we should use
            if ~isempty(delays) && ~isempty(dopplers) && ~isempty(gains)
                % take the assigned channels
                delays_len = length(delays);
                dopplers_len = length(dopplers);
                gains_len = length(gains);
                if ~isvector(delays) && ~isvector(dopplers) && ~isvector(gains) && (delays_len ~= dopplers_len || delays_len ~= gains_len)
                    error("The delays, dopplers and gains do not have the same length.");
                else
                    self.delay_taps = delays;
                    self.doppler_taps = dopplers;
                    self.chan_coef = gains;
                    self.taps_num = delays_len;
                end
            elseif p > 0 && lmax >= 1 && (kmax > 0 || kmax_frac > 0)
                % generate random channels
                % input check
                if p > (lmax + 1)*(2*kmax+1)
                    error("The path number must be less than lmax*(2*kmax+1) = %d", (lmax + 1)*(2*kmax+1));
                end
                if lmax >= self.nSubcarNum
                    error("The maximal delay index must be less than the subcarrier number.");
                end
                % input check - update kmax if using fractional Doppler
                if ~is_fractional_doppler
                    if kmax > floor(self.nTimeslotNum/2)
                        error("The maximal Doppler index must be less than the half of the timeslot number");
                    end
                else
                    if kmax_frac > floor(self.nTimeslotNum/2)
                        error("The maximal fractional Doppler index must be less than the half of the timeslot number");
                    end
                end
                lmin= 1;
                kmin = -kmax;
                taps_max = (kmax - kmin + 1)*(lmax - lmin + 1);
                % create delay options [lmin, lmin, lmin, lmin+1, lmin+1, lmin+1 ...]
                delay_taps_all = kron(lmin:lmax, ones(1, kmax - kmin + 1)); 
                % create Doppler options [kmin, kmin+1, kmin+2 ... kmax, kmin ...]
                Doppler_taps_all = repmat(kmin:kmax, 1, lmax - lmin + 1);
                % We select P paths from all possible paths; that is, we do the randperm(taps_max) and we choose the first P items
                taps_selected_idx = Channl_SelectRandomPathIdx(taps_max, p);
                % set channel info
                self.delay_taps = delay_taps_all(taps_selected_idx);
                % the 1st minimal delay is 0
                self.delay_taps(find(self.delay_taps == min(self.delay_taps), 1)) = 0;
                self.doppler_taps = Doppler_taps_all(taps_selected_idx);
                % add fractional Doppler
                if is_fractional_doppler
                    doppler_taps_k_max_pos_idx = self.doppler_taps == kmax;
                    doppler_taps_k_max_neg_idx = self.doppler_taps == -kmax;
                    doppler_taps_k_other_idx = abs(self.doppler_taps)~= kmax;
                    frac_range_max_pos = rand(1, p)*(kmax_frac - kmax + 0.5) - 0.5;
                    frac_range_max_neg = rand(1, p)*(kmax - kmax_frac - 0.5) + 0.5;
                    frac_range_others = rand(1, p) - 0.5;
                    frac_range_all = frac_range_max_pos.*doppler_taps_k_max_pos_idx + frac_range_max_neg.*doppler_taps_k_max_neg_idx + frac_range_others.*doppler_taps_k_other_idx;
                    self.doppler_taps = self.doppler_taps + frac_range_all;
                end
                self.chan_coef = sqrt(1/p)*(sqrt(1/2) * (randn(1, p)+1i*randn(1, p)));
                self.taps_num = p;
            else
                error("Channel Infomation is not recognised.");
            end
        end
        
        % add a path to the channel (this does not influence other existing paths)
        % @hi:      the path gain (linear gain)
        % @li:      the delay
        % @ki:      the Doppler shift
        function addChannelPath(self, hi, li, ki)
            if isempty(self.chan_coef)
                self.chan_coef = hi;
            else
                self.chan_coef = [self.chan_coef, hi];
            end
            if isempty(self.delay_taps)
                self.delay_taps = li;
            else
                self.delay_taps = [self.delay_taps, li];
            end
            if isempty(self.doppler_taps)
                self.doppler_taps = ki;
            else
                self.doppler_taps = [self.doppler_taps, ki];
            end
            self.taps_num = self.taps_num + 1;
        end
        
        % pass the channel
        % @No: noise power (a scalar) or a given noise vector
        function r = passChannel(self, No)
            % input check
            if isscalar(No)
                if No < 0
                    error("The noise power must be positive.");
                end
            elseif isvector(No)
                if length(No) ~= self.nSubcarNum*self.nTimeslotNum
                    error("The noise vector length must be %d.", self.nSubcarNum*self.nTimeslotNum);
                end
            else
                error("The noise input must be a scalar for power or a vector for fixed noise.");
            end
            
            % add CP
            cp_len = max(self.delay_taps);
            s_cp = Channel_AddCP(self.s, cp_len);
            % pass the channel
            s_chan = 0;
            for tap_id = 1:self.taps_num
                s_chan_tmp = self.chan_coef(tap_id)*circshift( ...
                    [ ...
                        s_cp.*exp(1j*2*pi/self.nSubcarNum*(-cp_len:-cp_len+length(s_cp)-1)*self.doppler_taps(tap_id)/self.nTimeslotNum).'; ...
                        zeros(cp_len,1) ...
                    ], ...
                    self.delay_taps(tap_id) ...
                );
                s_chan = s_chan+s_chan_tmp;
            end
            self.r = s_chan;
            % remove CP
            self.r = self.r(cp_len+1:cp_len+(self.nTimeslotNum*self.nSubcarNum));
            % add noise
            if isscalar(No)
                if No > 0
                    noise = sqrt(No/2)*(randn(self.nSubcarNum*self.nTimeslotNum, 1) + 1i*randn(self.nSubcarNum*self.nTimeslotNum, 1));
                    self.r = self.r + noise;
                end
            elseif isvector(No)
                self.r = self.r + No;
            end
            % return
            r = self.r;
        end
        
        %{
        % estimated the channel
        % @threshold: the threshold to detect a path
        %}
        function [gains, delays, dopplers] = estimateChannel(self, varargin)
            % optional inputs - register
            inPar = inputParser;
            addParameter(inPar, 'threshold', [], @(x) isempty(x)||isscalar(x)&&isnumeric(x));
            inPar.KeepUnmatched = true;     % Allow unmatched cases
            inPar.CaseSensitive = false;    % Allow capital or small characters
            parse(inPar, varargin{:});
            % optional inputs - assign
            threshold = inPar.Results.threshold;
            % input check
            % input check - insert pilots & guards
            if ~self.isInsertPilotsAndGuards()
                error("The pilot & guard allocation has to be done before channel estimation.");
            end
            % input check - demodulation
            if isempty(self.Y_DD) || isempty(self.y_DD)
                error("Demodulation has to be done before symbol detection.");
            end
            % input check - threshold
            if ~isempty(threshold) && threshold < 0
                error("The threshould must be non-negative.");
            end
            
            % estimate the channel
            pilot_len = length(self.pilots);
            if pilot_len == 1
                for delay_id = self.ce_ydd_delay_beg:self.ce_ydd_delay_end
                    for doppl_id = self.ce_ydd_doppl_beg:self.ce_ydd_doppl_end
                        pss_y = self.Y_DD(doppl_id, delay_id);
                        if abs(pss_y) > threshold
                            pss_beta = exp(1j*2*pi*(self.pilot_loc_delay_1st - 1)/self.nSubcarNum*(doppl_id - self.pilot_loc_doppl_1st)/self.nTimeslotNum);
                            self.ce_chan_coef(end+1) = pss_y/self.pilots/pss_beta;
                            self.ce_delay_taps(end+1) = delay_id - self.pilot_loc_delay_1st;
                            self.ce_doppler_taps(end+1) = doppl_id - self.pilot_loc_doppl_1st;
                        end
                    end
                end
            else
                warning("Multiple pilots CE is not supported.");
            end
            % return
            gains = self.ce_chan_coef;
            delays = self.ce_delay_taps;
            dopplers = self.ce_doppler_taps;
        end
        
        %{
        detect
        @detect_type:             detect type
        @csi_type:                detect CSI type
        @No:                      estimated noise power (linear)
        @constellation:           the constellation (a vector)
        @chan_coef:               the estimated channel gains
        @delay_taps:              the estimated channel delays
        @doppler_taps:            the estimated channel 
        @sym_map:                 false by default. If true, the output will be mapped to the constellation
        %}
        function symbols = detect(self, detect_type, csi_type, No, constellation, varargin)
            % optional inputs - register
            inPar = inputParser;
            addParameter(inPar,"chan_coef", [], @isnumeric);
            addParameter(inPar,"delay_taps", [], @isnumeric);
            addParameter(inPar,"doppler_taps", [], @isnumeric);
            addParameter(inPar,"sym_map", false, @(x) isscalar(x)&islogical(x)); 
            inPar.KeepUnmatched = true;
            inPar.CaseSensitive = false;
            parse(inPar, varargin{:});
            % optional inputs - assign
            chan_coef_est = inPar.Results.chan_coef;
            delay_taps_est = inPar.Results.delay_taps;
            doppler_taps_est = inPar.Results.doppler_taps;
            sym_map = inPar.Results.sym_map;
            % input check
            % input check - demodulation
            if isempty(self.Y_DD) || isempty(self.y_DD)
                error("Demodulation has to be done before symbol detection.");
            end
            % input check - detect_type
            if ~isscalar(detect_type) && ~ismember(detect_type, self.DETECT_TYPES)
                error("The detctor type is illegal.");
            end
            % input check - detect_csi_type
            if ~isscalar(csi_type) && ~ismember(csi_type, self.DETECT_CSI_TYPES)
                error("The channel state information source is illegal.");
            end
            if csi_type == self.DETECT_CSI_CE
                if ~self.isInsertPilotsAndGuards()
                    error("The pilot & guard allocation has to be done before channel estimation.");
                elseif isempty(self.ce_delay_taps) || isempty(self.ce_doppler_taps) || isempty(self.ce_chan_coef)
                    error("Channel estimation has to be done before symbol detection.");
                end
            elseif csi_type == self.DETECT_CSI_IN
                if ~isvector(chan_coef_est) || ~isvector(delay_taps_est) || ~isvector(doppler_taps_est)
                    error("The input CSI (gains, delays and dopplers) must be vectors.");
                end
                chan_coef_est_len = length(chan_coef_est);
                if chan_coef_est_len ~= length(delay_taps_est) || chan_coef_est_len ~= length(doppler_taps_est)
                    error("The input CSI (gains, delays and dopplers) must have the same length.");
                end
            end
            % input check - No
            if ~isscalar(No)
                error("The noise power(linear) must be a scalar.");
            end
            if ~isvector(constellation)
                error("The constellation must be a vector.");
            else
                % to row vector
                self.constellation = constellation(:).';
                self.constellation_len = length(constellation);
            end
            
            % retrieve the channel information
            if csi_type == self.DETECT_CSI_PERFECT
                chan_coef_est = self.chan_coef;
                delay_taps_est = self.delay_taps;
                doppler_taps_est = self.doppler_taps;
            end
            if csi_type ==  self.DETECT_CSI_CE
                chan_coef_est = self.ce_chan_coef;
                delay_taps_est = self.ce_delay_taps;
                doppler_taps_est = self.ce_doppler_taps;
            end
            % estimate the symbols
            switch detect_type
                case self.DETECT_MP_BASE
                    symbols = self.detectMPBase(No, chan_coef_est, delay_taps_est, doppler_taps_est);
            end
            % hard estimation
            if sym_map
                if detect_type ~= self.DETECT_MP_BASE
                    symbols = self.symmap(symbols);
                end
            end
        end
    end
    
    %% Getters & Setters
    methods
        %{
        Get the channel matrix in Delay Doppler Domain (using the rectangular waveform)
        @chan_coef:         the channel gains
        @delay_taps:        the channel delays
        @doppler_taps:      the channel Dopplers
        @only_for_data:     whether the channel is only for data (by default true). If you want to get the entire H_DD when using pilos and/or guards, you should manullay set it to false.
        %}
        function H_DD = getChannel(self, varargin)
            % optional inputs - register
            inPar = inputParser;
            addParameter(inPar,"chan_coef", [], @(x) isvector(x)&&isnumeric(x));
            addParameter(inPar,"delay_taps", [], @(x) isvector(x)&&isnumeric(x));
            addParameter(inPar,"doppler_taps", [], @(x) isvector(x)&&isnumeric(x));
            addParameter(inPar,"only_for_data", true, @(x) isscalar(x)&&islogical(x));
            inPar.KeepUnmatched = true;
            inPar.CaseSensitive = false;
            parse(inPar, varargin{:});
            % optional inputs - assign
            in_chan_coef = inPar.Results.chan_coef;
            in_delay_taps = inPar.Results.delay_taps;
            in_doppler_taps = inPar.Results.doppler_taps;
            only_for_data = inPar.Results.only_for_data;
            % input check
            % input check - csi
            csi_cond = sum(~isempty(in_chan_coef) + ~isempty(in_delay_taps) + ~isempty(in_doppler_taps));
            if csi_cond == 1 || csi_cond == 2
                error("If you give CSI, you have to give channel gains, delays and Dopplers together.");
            end
            in_taps_num = 0;
            if csi_cond == 3
                in_taps_num = length(in_chan_coef);
                if in_taps_num ~= length(in_delay_taps) && in_taps_num ~= length(in_doppler_taps)
                    error("The input CSI (gains, delays and dopplers) must have the same length.");
                end
            else
                in_taps_num = self.taps_num;
            end
            
            % reset the CSI if not input
            if csi_cond == 0
                in_chan_coef = self.chan_coef;
                in_delay_taps = self.delay_taps;
                in_doppler_taps = self.doppler_taps;
            end
            % build H_DD
            H_DD = zeros(self.nTimeslotNum*self.nSubcarNum, self.nTimeslotNum*self.nSubcarNum); % intialize the return channel
            dftmat = dftmtx(self.nTimeslotNum)/sqrt(self.nTimeslotNum); % DFT matrix
            idftmat = conj(dftmat); % IDFT matrix 
            piMat = eye(self.nTimeslotNum*self.nSubcarNum); % permutation matrix (from the delay) -> pi
            % accumulate all paths
            for tap_id = 1:in_taps_num
                hi = in_chan_coef(tap_id);
                li = in_delay_taps(tap_id);
                ki = in_doppler_taps(tap_id);
                % delay
                piMati = circshift(piMat, li); 
                % Doppler
                deltaMat_diag = exp(1j*2*pi*ki/(self.nTimeslotNum*self.nSubcarNum)*(0:1:self.nTimeslotNum*self.nSubcarNum-1));
                deltaMati = diag(deltaMat_diag);
                % Pi, Qi, & Ti
                Pi = kron(dftmat, eye(self.nSubcarNum))*piMati; 
                Qi = deltaMati*kron(idftmat, eye(self.nSubcarNum));
                Ti = Pi*Qi;
                H_DD = H_DD + hi*Ti;
            end
            % remove redundant values
            if self.isInsertPilotsAndGuards() && only_for_data
                invalud_row = NaN(1, self.nSubcarNum*self.nTimeslotNum);
                invalud_col = NaN(self.nSubcarNum*self.nTimeslotNum, 1);
                % mark redundant values - columns (X_DD invalid)
                for doppl_id = self.ce_xdd_doppl_beg:self.ce_xdd_doppl_end
                    for delay_id = self.ce_xdd_delay_beg:self.ce_xdd_delay_end
                        col_id = (doppl_id-1)*self.nSubcarNum + delay_id;
                        H_DD(:, col_id) = invalud_col;
                        if doppl_id == 2 && delay_id == 50
                            assert(sum(isnan(H_DD(:, col_id))) == self.nSubcarNum*self.nTimeslotNum);
                        end
                    end
                end
                % mark redundant values - rows (Y_DD invalid)
                for doppl_id = self.ce_ydd_doppl_beg:self.ce_ydd_doppl_end
                    for delay_id = self.ce_ydd_delay_beg:self.ce_ydd_delay_end
                        row_id = (doppl_id-1)*self.nSubcarNum + delay_id;
                        H_DD(row_id, :) = invalud_row;
                    end
                end
                % remove redundant values
                % remove redundant values - columns
                col_idx = ((sum(isnan(H_DD)) == self.nSubcarNum*self.nTimeslotNum));
                H_DD(:, col_idx) = [];
                % remove - rows
                [~, H_DD_col_num] = size(H_DD);
                row_idx = sum(isnan(H_DD), 2) == H_DD_col_num;
                H_DD(row_idx, :) = [];
            end
        end
        
        %{
        get the channel state information
        @sort_by_gain: sort axis
        @sort_by_delay_doppler: sort axes
        @sort_by_doppler_delay: sort axes
        @descend: sort direction
        %}
        function [gains, delays, dopplers] = getCSI(self, varargin)
            % optional inputs - register
            inPar = inputParser;
            addParameter(inPar,"sort_by_gain", false, @(x) isscalar(x)&&islogical(x));
            addParameter(inPar,"sort_by_delay_doppler", false, @(x) isscalar(x)&&islogical(x));
            addParameter(inPar,"sort_by_doppler_delay", false, @(x) isscalar(x)&&islogical(x));
            addParameter(inPar,"descend", false, @(x) isscalar(x)&&islogical(x));
            inPar.KeepUnmatched = true;
            inPar.CaseSensitive = false;
            parse(inPar, varargin{:});
            % optional inputs - assign
            sort_by_gain = inPar.Results.sort_by_gain;
            sort_by_delay_doppler = inPar.Results.sort_by_delay_doppler;
            sort_by_doppler_delay = inPar.Results.sort_by_doppler_delay;
            descend = inPar.Results.descend;
            % input check
            if sort_by_gain + sort_by_delay_doppler + sort_by_doppler_delay > 1
                error("Cannot sort following over two orders.");
            end
            
            % retrieve CSI
            gains = self.chan_coef(:).';
            delays = self.delay_taps(:).';
            dopplers = self.doppler_taps(:).';
            % only sort when there are multiple taps
            if self.taps_num > 1
                if sort_by_gain
                    [gains, sort_idx] = sort(gains, sort_direction);
                    delays = delays(sort_idx);
                    dopplers = dopplers(sort_idx);
                elseif sort_by_delay_doppler
                    arrs_sorted = self.sortArrs([gains;delays;dopplers], "first_order_arr", 2, "second_order_arr", 3, "descend", descend);
                    gains = arrs_sorted(1, :);
                    delays = arrs_sorted(2, :);
                    dopplers = arrs_sorted(3, :);
                elseif sort_by_doppler_delay
                    arrs_sorted = self.sortArrs([gains;delays;dopplers], "first_order_arr", 3, "second_order_arr", 2, "descend", descend);
                    gains = arrs_sorted(1, :);
                    delays = arrs_sorted(2, :);
                    dopplers = arrs_sorted(3, :);
                end
            end
        end
        
        
        %{
        get the signal in the TF domain
        %}
        function X_TF = getX2TF(self)
            X_TF = self.X_TF;
        end
        
        %{
        get the signal in the time domain
        %}
        function s = getX2T(self, varargin)
            % Inputs Name-Value Pair 
            inPar = inputParser;
            addParameter(inPar,'fft_size', 0, @isnumeric);
            inPar.KeepUnmatched = true;     % Allow unmatched cases
            inPar.CaseSensitive = false;    % Allow capital or small characters
            % fft size
            parse(inPar, varargin{:}); 
            fft_size = inPar.Results.fft_size;
            
            % if fft resolution is lower than subcarrier number, we choose the subcarrier number as the resolution
            if fft_size < self.nSubcarNum
                s = self.s;
            else 
                % Heisenberg transform
                s_mat = ifft(self.X_TF, fft_size)*sqrt(self.nSubcarNum);
                % vectorize
                s = s_mat(:);
            end
        end
        
        %{
        get the received signal in delay Doppler domain
        %}
        function Y_DD = getYDD(self)
            Y_DD = self.Y_DD;
        end
    end
    
    %% detectors
    methods(Access=private)
        %{
        MP base (proposed by P. Raviteja in 2017) from Emanuele Viterbo Research Group
        @No:              the estimated noise power
        @constellation:   the constellation (a vector)
        @chan_coef:       the channel coefficient
        @delay_taps:      the delay indices
        @Doppler_taps:    the Doppler indices
        @n_ite:           the iteration number (200 by default)
        @delta_fra:       the percentage for taking the values in the current iteration
        %}
        function symbols = detectMPBase(self, No, chan_coef, delay_taps, Doppler_taps, varargin)
            % input check
            if ~isvector(chan_coef)
                error("The channel coefficient must be a vector.");
            end
            if ~isvector(delay_taps)
                error("The delay indices must be a vector.");
            end
            if ~isvector(Doppler_taps)
                error("The Doppler indices must be a vector.");
            end
            taps = length(chan_coef);
            if taps ~= length(delay_taps) || taps ~= length(Doppler_taps)
                error("The channel parameters do not have the same length.");
            end
            % optional inputs - register
            default_n_ite = 200;
            default_delta_fra = 0.6;
            inPar = inputParser;
            addParameter(inPar,"n_ite", default_n_ite, @(x) isscalar(x)&&isnumeric(x));
            addParameter(inPar,"delta_fra", default_delta_fra, @(x) isscalar(x)&&isnumeric(x));
            inPar.KeepUnmatched = true;
            inPar.CaseSensitive = false;
            parse(inPar, varargin{:});
            % optional inputs
            n_ite = inPar.Results.n_ite;
            delta_fra = inPar.Results.delta_fra;
            % init
            yv = reshape(self.Y_DD, self.nTimeslotNum*self.nSubcarNum, 1);
            % init
            % init - observation nodes y[d]
            mean_int = zeros(self.nTimeslotNum*self.nSubcarNum,taps);
            var_int = zeros(self.nTimeslotNum*self.nSubcarNum,taps);
            % init - P_c,d (P_tap, d): y[d]<-x[c], aj, i.e., y[d]<-x[tap_i], aj.
            p_map = ones(self.nTimeslotNum*self.nSubcarNum,taps, self.constellation_len)*(1/self.constellation_len);
            % detect
            conv_rate_prev = -0.1;
            for ite=1:n_ite
                % Update mean and var (in the view of y[d])
                for ele1=1:1:self.nSubcarNum
                    for ele2=1:1:self.nTimeslotNum
                        % CE - jump the area for channel estimation
                        if self.isInsertPilotsAndGuards()
                            if ele1 >= self.ce_ydd_delay_beg && ele1 <= self.ce_ydd_delay_end && ele2 >= self.ce_ydd_doppl_beg && ele2 <= self.ce_ydd_doppl_end
                                continue;
                            end
                        end
                        % origianl part
                        mean_int_hat = zeros(taps,1);
                        var_int_hat = zeros(taps,1);
                        for tap_no=1:taps
                            m = ele1-1-delay_taps(tap_no)+1;
                            add_term = exp(1i*2*(pi/self.nSubcarNum)*(m-1)*(Doppler_taps(tap_no)/self.nTimeslotNum));
                            % calculate the α_i(k,l) in the paper
                            mp_ch_dd_alpha = 1;
                            if ele1-1<delay_taps(tap_no)
                                n = mod(ele2-1-Doppler_taps(tap_no),self.nTimeslotNum) + 1;
                                mp_ch_dd_alpha = exp(-1i*2*pi*((n-1)/self.nTimeslotNum));
                            end
                            new_chan = add_term * (mp_ch_dd_alpha) * chan_coef(tap_no);
                            % calculate the mean and variance from other taps contributing to this tap
                            for i2=1:1:self.constellation_len
                                % jump for CE area if use CE
                                if self.isInsertPilotsAndGuards()
                                    % only consider the place we have values
                                    mp_ob_mean_delay_sour = mod(ele1 + delay_taps(tap_no), self.nSubcarNum) + 1;
                                    mp_ob_mean_doppl_sour = mod(ele2 + Doppler_taps(tap_no), self.nTimeslotNum) + 1;
                                    if mp_ob_mean_delay_sour >= self.ce_xdd_delay_beg && mp_ob_mean_delay_sour <= self.ce_xdd_delay_end && mp_ob_mean_doppl_sour >= self.ce_xdd_doppl_beg && mp_ob_mean_doppl_sour <= self.ce_xdd_doppl_end
                                        continue;
                                    end
                                end
                                mean_int_hat(tap_no) = mean_int_hat(tap_no) + p_map(self.nTimeslotNum*(ele1-1)+ele2,tap_no,i2) * self.constellation(i2);
                                var_int_hat(tap_no) = var_int_hat(tap_no) + p_map(self.nTimeslotNum*(ele1-1)+ele2,tap_no,i2) * abs(self.constellation(i2))^2;
                            end
                            mean_int_hat(tap_no) = mean_int_hat(tap_no) * new_chan;
                            var_int_hat(tap_no) = var_int_hat(tap_no) * abs(new_chan)^2;
                            var_int_hat(tap_no) = var_int_hat(tap_no) - abs(mean_int_hat(tap_no))^2;
                        end

                        mean_int_sum = sum(mean_int_hat);
                        var_int_sum = sum(var_int_hat)+(No);
                        % remove the c tap from Gaussian variable estimation sum, so it remains the 
                        for tap_no=1:taps
                            mean_int(self.nTimeslotNum*(ele1-1)+ele2,tap_no) = mean_int_sum - mean_int_hat(tap_no);
                            var_int(self.nTimeslotNum*(ele1-1)+ele2,tap_no) = var_int_sum - var_int_hat(tap_no);
                        end

                    end
                end
                % Update probabilities (in the view of x[c])
                sum_prob_comp = zeros(self.nTimeslotNum*self.nSubcarNum, self.constellation_len);
                dum_eff_ele1 = zeros(taps,1);
                dum_eff_ele2 = zeros(taps,1);
                for ele1=1:1:self.nSubcarNum
                    for ele2=1:1:self.nTimeslotNum
                        % CE - jump the area for channel estimation
                        if self.isInsertPilotsAndGuards()
                            if ele1 >= self.ce_xdd_delay_beg && ele1 <= self.ce_xdd_delay_end && ele2 >= self.ce_xdd_doppl_beg && ele2 <= self.ce_xdd_doppl_end
                                continue;
                            end
                        end
                        
                        % original code
                        dum_sum_prob = zeros(self.constellation_len,1);
                        log_te_var = zeros(taps,self.constellation_len);
                        for tap_no=1:taps
                            if ele1+delay_taps(tap_no)<=self.nSubcarNum
                                eff_ele1 = ele1 + delay_taps(tap_no);
                                add_term = exp(1i*2*(pi/self.nSubcarNum)*(ele1-1)*(Doppler_taps(tap_no)/self.nTimeslotNum));
                                int_flag = 0;
                            else
                                eff_ele1 = ele1 + delay_taps(tap_no)- self.nSubcarNum;
                                add_term = exp(1i*2*(pi/self.nSubcarNum)*(ele1-1-self.nSubcarNum)*(Doppler_taps(tap_no)/self.nTimeslotNum));
                                int_flag = 1;
                            end
                            add_term1 = 1;
                            if int_flag==1
                                add_term1 = exp(-1i*2*pi*((ele2-1)/self.nTimeslotNum));
                            end
                            eff_ele2 = mod(ele2-1+Doppler_taps(tap_no),self.nTimeslotNum) + 1;
                            dum_eff_ele1(tap_no) = eff_ele1;
                            dum_eff_ele2(tap_no) = eff_ele2;
                            new_chan = add_term * add_term1 * chan_coef(tap_no);
                            
                            % only consider the place we have values
                            mp_va_mean_delay_sour = mod(ele1 + delay_taps(tap_no), self.nSubcarNum);
                            mp_va_mean_doppl_sour = mod(ele2 + Doppler_taps(tap_no), self.nTimeslotNum) + 1;
                            if mp_va_mean_delay_sour >= self.ce_xdd_delay_beg && mp_va_mean_delay_sour <= self.ce_xdd_delay_end && mp_va_mean_doppl_sour >= self.ce_xdd_doppl_beg && mp_va_mean_doppl_sour <= self.ce_xdd_doppl_end
                                continue;
                            end
                              
                            for i2=1:1:self.constellation_len
                                dum_sum_prob(i2) = abs(yv(self.nTimeslotNum*(eff_ele1-1)+eff_ele2)- mean_int(self.nTimeslotNum*(eff_ele1-1)+eff_ele2,tap_no) - new_chan * self.constellation(i2))^2;
                                dum_sum_prob(i2)= -(dum_sum_prob(i2)/var_int(self.nTimeslotNum*(eff_ele1-1)+eff_ele2,tap_no));
                            end
                            dum_sum = dum_sum_prob - max(dum_sum_prob);
                            dum1 = sum(exp(dum_sum));
                            log_te_var(tap_no,:) = dum_sum - log(dum1);
                        end
                        for i2=1:1:self.constellation_len
                            ln_qi(i2) = sum(log_te_var(:,i2));
                        end
                        dum_sum = exp(ln_qi - max(ln_qi));
                        dum1 = sum(dum_sum);
                        sum_prob_comp(self.nTimeslotNum*(ele1-1)+ele2,:) = dum_sum/dum1;
                        for tap_no=1:1:taps
                            eff_ele1 = dum_eff_ele1(tap_no);
                            eff_ele2 = dum_eff_ele2(tap_no);
                            
                            if eff_ele1 >= self.ce_xdd_delay_beg && eff_ele1 <= self.ce_xdd_delay_end && eff_ele2 >= self.ce_xdd_doppl_beg && eff_ele2 <= self.ce_xdd_doppl_end
                                continue;
                            end

                            dum_sum = log_te_var(tap_no,:);
                            ln_qi_loc = ln_qi - dum_sum;
                            dum_sum = exp(ln_qi_loc - max(ln_qi_loc));
                            dum1 = sum(dum_sum);
                            p_map(self.nTimeslotNum*(eff_ele1-1)+eff_ele2,tap_no,:) = (dum_sum/dum1)*delta_fra + (1-delta_fra)*reshape(p_map(self.nTimeslotNum*(eff_ele1-1)+eff_ele2,tap_no,:),1,self.constellation_len);
                        end
                    end
                end
                % convergence indicator
                conv_rate =  sum(max(sum_prob_comp,[],2)>0.99)/(self.nTimeslotNum*self.nSubcarNum);
                if conv_rate==1
                    sum_prob_fin = sum_prob_comp;
                    break;
                elseif conv_rate > conv_rate_prev
                    conv_rate_prev = conv_rate;
                    sum_prob_fin = sum_prob_comp;
                elseif (conv_rate < conv_rate_prev - 0.2) && conv_rate_prev > 0.95
                    break;
                end
            end
            % estimate symbols
            X_DD_est = zeros(self.nTimeslotNum, self.nSubcarNum);
            for ele1=1:1:self.nSubcarNum
                for ele2=1:1:self.nTimeslotNum
                    % CE - jump the area for channel estimation
                    if self.isInsertPilotsAndGuards()
                        if ele1 >= self.ce_xdd_delay_beg && ele1 <= self.ce_xdd_delay_end && ele2 >= self.ce_xdd_doppl_beg && ele2 <= self.ce_xdd_doppl_end
                            continue;
                        end
                    end
                    [~,pos] = max(sum_prob_fin(self.nTimeslotNum*(ele1-1)+ele2,:));
                    X_DD_est(ele2,ele1) = self.constellation(pos);
                end
            end
            % extract symbols
            symbols = zeros(self.nTimeslotNum*self.nSubcarNum - self.ce_xdd_num, 1);
            % return
            if ~self.isInsertPilotsAndGuards()
                % no CE
                symbols = X_DD_est.';
                symbols = symbols(:);
            else
                % CE in use 
                symbols_id = 1;
                for doppl_id = 1:self.nTimeslotNum
                    for delay_id = 1:self.nSubcarNum
                        if doppl_id<self.ce_xdd_doppl_beg || doppl_id>self.ce_xdd_doppl_end || delay_id<self.ce_xdd_delay_beg || delay_id>self.ce_xdd_delay_end
                            symbols(symbols_id) = X_DD_est(doppl_id, delay_id);
                            symbols_id = symbols_id + 1;
                        end
                    end
                end
            end
        end
    end
    
    %% support functions - general
    methods(Access=private)
        %{
        check whether the pilots & guards are inserted or not
        %}
        function is_done = isInsertPilotsAndGuards(self)
            % check whether pilots is assigned or not
            is_done = ~isempty(self.pilots);
            % check whether X_DD_invalid area is calculated
            is_done = is_done && ~isnan(self.ce_xdd_delay_beg) && ~isnan(self.ce_xdd_delay_end) && ~isnan(self.ce_xdd_doppl_beg) && ~isnan(self.ce_xdd_doppl_end);
            % check whether CE area is calulated
            is_done = is_done && ~isnan(self.ce_ydd_delay_beg) && ~isnan(self.ce_ydd_delay_end) && ~isnan(self.ce_ydd_doppl_beg) && ~isnan(self.ce_ydd_doppl_end);
        end

        %{
        check whether the current position is in channel estimation area
        @ce_area_tag:   notify whether CE area to check
        @pos_delay:     the position on the delay axis for matrix or th position on the Doppler-delay axis for the vector
        @pos_doppl:     the position on the Doppler axis. Not given means the position is for a Doppler-delay vector
        %}
        function is_in = isCurPosInCEXDDArea(delay_pos, varargin)
            is_in = self.isCurPosInCEArea(self.CE_AREA_TAG_XDD, delay_pos, varargin);
        end
        function is_in = isCurPosInCEYDDArea(delay_pos, varargin)
            is_in = self.isCurPosInCEArea(self.CE_AREA_TAG_YDD, delay_pos, varargin);
        end
        function is_in = isCurPosInCEArea(self, ce_area_tag, delay_pos, varargin)
            % input check
            if ismember(ce_area_tag, self.CE_AREA_TAGS)
                error("The CE area tag is illegal.");
            end
            if isempty(varargin)
                if delay_pos <= 0 || delay_pos > self.nSubcarNum*self.nTimeslotNum
                    error("The vector position is out of the OTFS size.");
                end
            else
                if delay_pos <= 0 || delay_pos > self.nSubcarNum
                    error("The delay position is out of the subcarrier number.");
                end
                if varargin{1} <= 0 || varargin{1} > self.nTimeslotNum
                    error("The Doppler position is out of the timeslot number.");
                end
            end
            
            % recalculate the position
            pos_doppl = 0;
            if isempty(varargin)
                pos_doppl = delay_pos/self.nSubcarNum;
                if pos_doppl == floor(pos_doppl)
                    pos_doppl = floor(pos_doppl);
                else
                    pos_doppl = floor(pos_doppl) + 1;
                end
                delay_pos = delay_pos - (pos_doppl-1)*self.nSubcarNum;
            else
                pos_doppl = varargin{1};
            end
            % decide
            if ce_area_tag == self.CE_AREA_TAG_XDD
                is_in = delay_pos >= self.ce_xdd_delay_beg && delay_pos <= self.ce_xdd_delay_end && pos_doppl >= self.ce_xdd_doppl_beg && pos_doppl <= self.ce_xdd_doppl_end;
            end
            if ce_area_tag == self.CE_AREA_TAG_YDD
                is_in = delay_pos >= self.ce_ydd_delay_beg && delay_pos <= self.ce_ydd_delay_end && pos_doppl >= self.ce_ydd_doppl_beg && pos_doppl <= self.ce_ydd_doppl_end;
            end
        end

        %{
        symbol mapping (hard)
        @syms: a vector of symbols
        @constellation: the constellation to map
        %}
        function syms_mapped = symmap(self, syms)
            if ~isvector(syms)
                error("Symbols must be into a vector form to map.");
            end
            % the input must be a column vector
            is_syms_col = iscolumn(syms);
            syms = syms(:);
            syms_dis = abs(syms - self.constellation).^2;
            [~,syms_dis_min_idx] =  min(syms_dis,[],2);
            syms_mapped = self.constellation(syms_dis_min_idx);
            if is_syms_col
                syms_mapped = syms_mapped(:);
            end
        end
    end
        
    %% support function - Maths+
    methods(Access=private)
        %{
        sort arrays based on the 
        @arrs:                arrays [arr_num, arr_data_num]
        @first_order_arr:     the 1st array id we follow the order
        @second_order_arr:    the 2nd array id we follow the order
        @descend:             descend order
        %}
        function arrs = sortArrs(self, arrs, varargin)
            % optional inputs - register
            inPar = inputParser;
            addParameter(inPar,"first_order_arr", false, @(x) isscalar(x)&&isnumeric(x));
            addParameter(inPar,"second_order_arr", false, @(x) isscalar(x)&&isnumeric(x));
            addParameter(inPar,"descend", false, @(x) isscalar(x)&&islogical(x));
            inPar.KeepUnmatched = true;
            inPar.CaseSensitive = false;
            parse(inPar, varargin{:});
            % optional inputs - assign
            first_order_arr = inPar.Results.first_order_arr;
            second_order_arr = inPar.Results.second_order_arr;
            descend = inPar.Results.descend;
            % input check
            % input check - arrays
            if isvector(arrs) || ~ismatrix(arrs)
                error("Arrays must be given in a matrix form [arr_num, arr_data_num].");
            end
            [arr_num, arr_data_num] = size(arrs);
            % input check - first_order_arr & second_order_arr
            if first_order_arr < 0 || first_order_arr > arr_num || second_order_arr < 0 || second_order_arr > arr_num || first_order_arr == second_order_arr
                error("The selected two arrays are illegal.");
            end
            
            % retrieve 1st and 2nd arr
            arr1 = arrs(first_order_arr, :);
            arr2 = arrs(second_order_arr, :);
            % sort - order
            sort_direction = 'ascend';
            if descend
                sort_direction = 'descend';
            end
            % sort - 1st arr
            [~, sort_idx_dim1] = sort(arr1, sort_direction);
            for arr_id = 1:arr_num
                cur_arr = arrs(arr_id, :);
                arrs(arr_id, :) = cur_arr(sort_idx_dim1);
            end
            % sort - 2nd arr
            arr1 = arrs(first_order_arr, :);
            arr2 = arrs(second_order_arr, :);
            arr_data_id = 1;
            while arr_data_id <= arr_data_num
                % find the beginning and the end of the segment
                seg_beg = find(arr1 == arr1(arr_data_id), 1, "first");
                seg_end = find(arr1 == arr1(arr_data_id), 1, "last");
                % move to the next segment
                arr_data_id = seg_end + 1;
                % jump if there is only 1 element
                if seg_beg == seg_end
                    continue;
                end
                % segment - retrieve
                seg_arr2 = arr2(seg_beg:seg_end);
                % segment - sort
                [~, sort_idx_dim2] = sort(seg_arr2, sort_direction);
                % segments - sort
                segs = arrs(:, seg_beg:seg_end);
                for arr_id = 1:arr_num
                    segs_row = segs(arr_id, :);
                    segs(arr_id, :) = segs_row(sort_idx_dim2);
                end
                % arrs - sort
                arrs(:, seg_beg:seg_end) = segs;
            end
        end
    end
end

%% independent functions
% select p random path indices from given paths
% @taps_max: the maximal path number
% @p: the path number we select
function taps_selected_idx = Channl_SelectRandomPathIdx(taps_max, p)
    taps_idx_chaotic = randperm(taps_max);
    taps_selected_idx = taps_idx_chaotic(1:p);
end

% add cyclic prefix
% @s: time domain signal
% @cp_len: the cyclic prefix length
function s_cp = Channel_AddCP(s, cp_len)
    syms_len = length(s);
    s_cp = [s(syms_len - cp_len + 1 : syms_len);s];
end
