classdef OTFS < handle
    % constants
    properties(Constant)
        % Pulse Types
        PULSE_IDEAL = 10;   % ideal pulses (if we use ideal pulses, `ISI_CANCEL_CP` is forced to chosen)
        PULSE_RECTA = 20;   % rectangular pulses
        % cyclic prefix
        CP_ZERO = 10;                                % use zero guards to avoid ISI
        CP_ONE_FRAM = 20;                            % one cp for entire OTFS frame
        CP_ONE_FRAM_SUB = 21;                        % one cp for each OTFS subframe
        % Detect
        DETECT_MP_BASE = 10;                        % base OTFS MP detector proposed by P. Raviteja in 2018
        DETECT_TYPES = [OTFS.DETECT_MP_BASE];
        % CSI
        DETECT_CSI_PERFECT = 10;                    % perfect CSI
        DETECT_CSI_CE = 20;                         % CSI from channel estimation (its type is based on pilot type)
        DETECT_CSI_IN = 30;                         % CSI from the input
        DETECT_CSI_TYPES = [OTFS.DETECT_CSI_PERFECT, OTFS.DETECT_CSI_CE, OTFS.DETECT_CSI_IN];
    end
    properties
        % RG infomation
        nSubcarNum {mustBeInteger}                  % subcarrier number
        nTimeslotNum {mustBeInteger}                % timeslot number
        sig_len = 0;                                % the total signal length
        rg
        % OTFS process signal
        X_TF                                        % Tx value in the time-frequency(TF) domain
        s                                           % Tx value in the time domain (array)
        H                                           % channel in the time domain
        r                                           % Rx value in the time domain (array)
        Y_TF                                        % Rx value in the TF domain
        Y_DD
        % channel
        taps_num = 0                                % paths number           
        chan_coef                                   % path gain, a row vector
        delay_taps                                  % delay index, a row vector
        doppler_taps                                % doppler index (integers or fractional numbers), a row vector
        % pulse
        pulse_type = OTFS.PULSE_RECTA;
        % cp
        cp_type = OTFS.CP_ONE_FRAM;
        cp_len = 0;
        % detection
        constellation                               % constellation
        constellation_len = 0;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % General OTFS Methods
    methods
        %{
        modulate (use fast method by default)
        @rg:        an OTFS resource grid
        @isFast:    DD domain -> TD domain (no X_TF) 
        %}
        function modulate(self, rg, varargin)
            % optional inputs - register
            inPar = inputParser;
            addParameter(inPar,"isFast", true, @(x) isscalar(x)&islogical(x)); 
            inPar.KeepUnmatched = true;
            inPar.CaseSensitive = false;
            parse(inPar, varargin{:});
            isFast = inPar.Results.isFast;
            % input check
            if ~isa(rg, 'OTFSResGrid')
                error("The input must be an OTFS resource grid.");
            end
            % load RG
            [self.nSubcarNum, self.nTimeslotNum] = rg.getContentSize();
            self.sig_len = self.nSubcarNum*self.nTimeslotNum;
            if rg.isPulseIdeal()
                self.pulse_type = self.PULSE_IDEAL;
            elseif rg.isPulseRecta()
                self.pulse_type = self.PULSE_RECTA;
            else
                error("Pulse shaping is not given in the resource grid.");
            end
            self.rg = rg.clone();
            X_DD = self.rg.getContent();
            % modulate
            if self.pulse_type == self.PULSE_RECTA
                if isFast
                    s_mat = ifft(X_DD).'*sqrt(self.nTimeslotNum);
                    self.s = s_mat(:);
                else
                    X_FT = fft(ifft(X_DD).').'/sqrt(self.nSubcarNum/self.nTimeslotNum); % ISFFT 
                    self.X_TF = X_FT.'; % X_TF is [nSubcarNum, nTimeslotNum]
                    s_mat = ifft(self.X_TF)*sqrt(self.nSubcarNum); % Heisenberg transform
                    self.s = s_mat(:);
                end
            end
        end

        % set channel (in1, in2, in3)
        % set a fixed chanel (at least two paths, if you want to add one fixed path, call `setChannelExtra`)
        % @in1->his:        the path gains
        % @in2->lis:        the delays
        % @in3->kis:        the doppler shifts
        % set a random channel (overwritten the channel setting; use Rayleigh fading if not select channel model)
        % @in1->p:          the path number
        % @in2->lmax:       the maxmimal delay index
        % @in3->kmax:       the maximal Doppler index (can be fractional)
        % @force_frac:      use fractional Doppler (force)
        % @isAWGN:          use awgn
        % @isRician:        use Rician fading
        function setChannel(self, in1, in2, in3, varargin)
            % optional inputs - register
            inPar = inputParser;
            addParameter(inPar,"force_frac", false, @(x) isscalar(x)&islogical(x)); 
            addParameter(inPar,"isAWGN", false, @(x) isscalar(x)&islogical(x)); 
            addParameter(inPar,"isRician", true, @(x) isscalar(x)&islogical(x)); 
            inPar.KeepUnmatched = true;
            inPar.CaseSensitive = false;
            parse(inPar, varargin{:});
            force_frac = inPar.Results.force_frac;
            isAWGN = inPar.Results.isAWGN;
            isRician = inPar.Results.isRician;
            % set the channel
            if ~isscalar(in1) && ~isscalar(in2) && ~isscalar(in3)
                % no scalar, giving fixed path
                in1_len = length(in1);
                in2_len = length(in2);
                in3_len = length(in3);
                if ~isvector(in1) || ~isvector(in2) || ~isvector(in3)
                    error("The path gains, delays and dopplers must be 1D vectors.");
                elseif in1_len ~= in2_len || in1_len ~= in3_len
                    error("The path gains, delays and dopplers do not have the same length.");
                else
                    self.taps_num = in1_len;
                    self.chan_coef = in1;
                    self.delay_taps = in2;
                    self.doppler_taps = in3;
                    self.cp_len = max(self.delay_taps);
                end
            elseif isscalar(in1) && isscalar(in2) && isscalar(in3)
                % scalar, random paths
                p = in1;
                lmax = in2;
                kmax = in3;
                kmax_frac = kmax - floor(kmax);
                kmax = floor(kmax);
                if kmax_frac > eps
                    force_frac = true;
                else
                    if ~isempty(varargin)
                        force_frac = varargin{1};
                        if force_frac
                            kmax_frac = 0.5;
                        end
                    end
                end
                % input check
                if p > (lmax + 1)*(2*kmax+1)
                    error("The path number must be less than lmax*(2*kmax+1) = %d", (lmax + 1)*(2*kmax+1));
                end
                lmin= 1;
                kmin = -kmax;
                taps_max = (kmax - kmin + 1)*(lmax - lmin + 1);
                delay_taps_all = kron(lmin:lmax, ones(1, kmax - kmin + 1)); % create delay options [lmin, lmin, lmin, lmin+1, lmin+1, lmin+1 ...]
                Doppler_taps_all = repmat(kmin:kmax, 1, lmax - lmin + 1); % create Doppler options [kmin, kmin+1, kmin+2 ... kmax, kmin ...]
                taps_selected_idx = self.shufSelectTopNIdx(taps_max, p); % select P paths from all possible paths
                % CSI - delay
                self.delay_taps = delay_taps_all(taps_selected_idx);
                self.delay_taps(find(self.delay_taps == min(self.delay_taps), 1)) = 0;
                % CSI - doppler
                self.doppler_taps = Doppler_taps_all(taps_selected_idx);
                if force_frac   % add fractional Doppler
                    doppler_taps_k_max_pos_idx = self.doppler_taps == kmax;
                    doppler_taps_k_max_neg_idx = self.doppler_taps == -kmax;
                    doppler_taps_k_other_idx = abs(self.doppler_taps)~= kmax;
                    frac_range_max_pos = rand(1, p)*(kmax_frac - kmax + 0.5) - 0.5;
                    frac_range_max_neg = rand(1, p)*(kmax - kmax_frac - 0.5) + 0.5;
                    frac_range_others = rand(1, p) - 0.5;
                    frac_range_all = frac_range_max_pos.*doppler_taps_k_max_pos_idx + frac_range_max_neg.*doppler_taps_k_max_neg_idx + frac_range_others.*doppler_taps_k_other_idx;
                    self.doppler_taps = self.doppler_taps + frac_range_all;
                end
                % CSI - others
                self.taps_num = p;
                if isAWGN
                    self.chan_coef = sqrt(1/p)*(sqrt(1/2) * (ones(1, p)+1i*ones(1, p)));
                elseif isRician
                    self.chan_coef = sqrt(1/p)*(sqrt(1/2) * (ones(1, p)+1i*ones(1, p))) + sqrt(1/p)*(sqrt(1/2) * (randn(1, p)+1i*randn(1, p)));
                else
                    self.chan_coef = sqrt(1/p)*(sqrt(1/2) * (randn(1, p)+1i*randn(1, p))); % use Rayleigh fading by default
                end
                self.cp_len = max(self.delay_taps);
            else
                error("The given CSI is not recognised.");
            end
        end
        
        %{
        add a path to the channel (this does not influence other existing paths)
        @hi:      the path gain (linear gain)
        @li:      the delay
        @ki:      the Doppler shift
        %}
        function setChannelExtra(self, hi, li, ki)
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
            self.cp_len = max(self.delay_taps);
        end
        
        %{
        pass the channel
        @No: noise power (a scalar) or a given noise vector
        %}
        function r = passChannel(self, No)
            % input check
            if isscalar(No)
                if No < 0
                    error("The noise power must be positive.");
                end
            elseif isvector(No)
                if length(No) ~= self.sig_len
                    error("The noise vector length must be %d.", self.sig_len);
                end
            else
                error("The noise input must be a scalar for power or a vector for fixed noise.");
            end
            % pass based on the pulse type
            s_chan = 0;
            if self.pulse_type == self.PULSE_IDEAL
                H_DD = self.buildIdealChannel(self.taps_num, self.chan_coef, self.delay_taps, self.doppler_taps);
                s_chan = H_DD*self.rg.getContent("isVector", true);
            elseif self.pulse_type == self.PULSE_RECTA
                [s_cp, s_cp_t] = self.addCP();
                for tap_id = 1:self.taps_num
                    s_chan_tmp = self.chan_coef(tap_id)*circshift([ ...
                            s_cp.*exp(1j*2*pi*self.doppler_taps(tap_id)*s_cp_t/self.nSubcarNum/self.nTimeslotNum); ...
                            zeros(self.cp_len,1) ...
                        ], ...
                        self.delay_taps(tap_id) ...
                    );
                    s_chan = s_chan+s_chan_tmp;
                end
                s_chan = self.removeCP(s_chan);
            end
            % add noise
            if isscalar(No)
                noise = sqrt(No/2)*(randn(self.sig_len, 1) + 1i*randn(self.sig_len, 1));
                s_chan = s_chan + noise;
            elseif isvector(No)
                s_chan = s_chan + No;
            end
            % assign to Rx
            if self.pulse_type == self.PULSE_IDEAL
                self.Y_DD = reshape(s_chan, self.nSubcarNum, self.nTimeslotNum).';
            elseif self.pulse_type == self.PULSE_RECTA
                self.r = s_chan;    % to time domain
                r = s_chan;
            end
        end

        %{
        demodulate (use fast method by default)
        @isFast:    TD domain -> DD domain (no Y_TF) 
        %}
        function rg = demodulate(self, varargin)
            % optional inputs - register
            inPar = inputParser;
            addParameter(inPar,"isFast", true, @(x) isscalar(x)&islogical(x)); 
            inPar.KeepUnmatched = true;
            inPar.CaseSensitive = false;
            parse(inPar, varargin{:});
            isFast = inPar.Results.isFast;
            % demodulate
            if self.pulse_type == self.PULSE_RECTA
                r_mat = reshape(self.r, self.nSubcarNum, self.nTimeslotNum);
                if isFast
                    self.Y_DD = fft(r_mat.')/sqrt(self.nTimeslotNum);
                else 
                    self.Y_TF = fft(r_mat)/sqrt(self.nSubcarNum);  % Wigner transform (Y_TF in [nSubcarNum, nTimeslotNum])
                    Y_FT = self.Y_TF.';
                    self.Y_DD = ifft(fft(Y_FT).').'/sqrt(self.nTimeslotNum/self.nSubcarNum); % SFFT (Y_DD in [Doppler, delay] or [nTimeslotNum ,nSubcarNum])
                end
            end
            self.rg.setContent(self.Y_DD);
            rg = self.rg;
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
                if ~isempty(self.ce_delay_taps) || isempty(self.ce_doppler_taps) || isempty(self.ce_chan_coef)
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Getters & Setters
    methods
        %{
        Get the channel matrix in Delay Doppler Domain
        @his:   the channel gains
        @lis:   the channel delays
        @kis:   the channel Dopplers
        @data_only: whether the channel is only for data (by default true). If you want to get the entire H_DD when using pilos and/or guards, you should manullay set it to false.
        %}
        function H_DD = getChannel(self, varargin)
            % input check & init
            p = self.taps_num;
            his = self.chan_coef;
            lis = self.delay_taps;
            kis = self.doppler_taps;
            data_only = true;
            if ~isempty(varargin)
                % load optional inputs
                if length(varargin) >= 3
                    his = varargin{1};
                    lis = varargin{2};
                    kis = varargin{3};
                    p = length(his);
                    if p ~= length(lis) && p ~= length(kis)
                        error("The input CSI (gains, delays and dopplers) must have the same length.");
                    end
                end
                % load paired values
                inPar = inputParser;
                addParameter(inPar,"data_only", true, @(x) isscalar(x)&&islogical(x));
                inPar.KeepUnmatched = true;
                inPar.CaseSensitive = false;
                if length(varargin) <= 2
                    parse(inPar, varargin{:});
                else
                    parse(inPar, varargin{4:end});
                end
                data_only = inPar.Results.data_only;
            end
            % build the channel
            if self.pulse_type == self.PULSE_IDEAL
                H_DD = self.buildIdealChannel(p, his, lis, kis);
            elseif self.pulse_type == self.PULSE_RECTA
                H_DD = self.buildRectaChannel(p, his, lis, kis);
            end
            % remove the channel for PG & CE
            if data_only
                H_DD = self.removeNoDAChannel(H_DD);
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
        function X_TF = getXTF(self)
            X_TF = self.X_TF;
        end
        
        %{
        get the signal in the time domain
        @fft_size: the size of fft
        %}
        function s = getXT(self, varargin)
            % set fft size
            fft_size = 0;
            if ~isempty(varargin)
                fft_size = floor(varargin{1});
            end
            % if fft resolution is lower than subcarrier number, we choose the subcarrier number as the resolution
            if fft_size < self.nSubcarNum
                s = self.s;
            else 
                s_mat = ifft(self.X_TF, fft_size)*sqrt(self.nSubcarNum); % Heisenberg transform
                s = s_mat(:); % vectorize
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % detectors
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
                    [~,pos] = max(sum_prob_fin(self.nTimeslotNum*(ele1-1)+ele2,:));
                    X_DD_est(ele2,ele1) = self.constellation(pos);
                end
            end
            % no CE
            symbols = X_DD_est.';
            symbols = symbols(:);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % private methods
    methods(Access=private)
        %{
        shuffle and select top n elements' indices
        @total: 
        %}
        function idx = shufSelectTopNIdx(~, taps_max, p)
            taps_idx_chaotic = randperm(taps_max);
            idx = taps_idx_chaotic(1:p);
        end

        %{
        add cyclic prefix
        %}
        function [s_cp, s_cp_t] = addCP(self)
            if self.cp_type == self.CP_ZERO
                s_cp = self.s;
                s_cp_t = (0:self.nSubcarNum*self.nTimeslotNum-1)';
            elseif self.cp_type == self.CP_ONE_FRAM
                s_cp = [self.s(end - self.cp_len + 1 : end);self.s];
                s_cp_t = (-self.cp_len:self.nSubcarNum*self.nTimeslotNum-1)';
            elseif self.cp_type == self.CP_ONE_FRAM_SUB
                s_mat = reshape(self.s, self.nSubcarNum, self.nTimeslotNum);
                s_mat = [s_mat(end - self.cp_len + 1 : end, :); s_mat];
                s_cp = s_mat(:);
                s_cp_t_mat = repmat((-self.cp_len:self.nSubcarNum-1)', 1, self.nTimeslotNum);
                s_cp_t = s_cp_t_mat(:);
            end
        end

        %{
        remove cyclic prefix
        @s_chan: channel output
        %}
        function s_chan_rm = removeCP(self, s_chan)
            if self.cp_type == self.CP_ZERO
                s_chan_rm = s_chan;
            elseif self.cp_type == self.CP_ONE_FRAM
                s_chan_rm = s_chan(self.cp_len+1 : self.cp_len+self.nTimeslotNum*self.nSubcarNum);
            elseif self.cp_type == self.CP_ONE_FRAM_SUB
                s_chan_mat = reshape(s_chan, self.nSubcarNum, self.nTimeslotNum);
                s_chan_rm_mat = s_chan_mat(self.cp_len+1:self.cp_len+self.nTimeslotNum*self.nSubcarNum, :);
                s_chan_rm = s_chan_rm_mat(:);
            end
        end

        %{
        sort arrays based on the 
        @arrs:                arrays [arr_num, arr_data_num]
        @first_order_arr:     the 1st array id we follow the order
        @second_order_arr:    the 2nd array id we follow the order
        @descend:             descend order
        %}
        function arrs = sortArrs(~, arrs, varargin)
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

        %{
        remove non-data from DD channel
        @H_DD: the DD channel from
        %}
        function H_DD = removeNoDAChannel(self, H_DD)
            invalud_row = NaN(1, self.sig_len);
            invalud_col = NaN(self.sig_len, 1);
            [pg_num, pg_delay_beg, pg_delay_end, pg_doppl_beg, pg_doppl_end] = self.rg.getAreaPG();
            [ce_num, ce_delay_beg, ce_delay_end, ce_doppl_beg, ce_doppl_end] = self.rg.getAreaCE();
            % mark redundant values - columns (PG area)
            if pg_num > 0
                for doppl_id = pg_doppl_beg:pg_doppl_end
                    for delay_id = pg_delay_beg:pg_delay_end
                        col_id = (doppl_id-1)*self.nSubcarNum + delay_id;
                        H_DD(:, col_id) = invalud_col;
                    end
                end
            end
            % mark redundant values - rows (CE area)
            if ce_num > 0
                for doppl_id = ce_doppl_beg:ce_doppl_end
                    for delay_id = ce_delay_beg:ce_delay_end
                        row_id = (doppl_id-1)*self.nSubcarNum + delay_id;
                        H_DD(row_id, :) = invalud_row;
                    end
                end
            end
            % remove
            % remove - columns
            if pg_num > 0
                col_idx = ((sum(isnan(H_DD)) == self.sig_len));
                H_DD(:, col_idx) = [];
            end
            % remove - rows
            if ce_num > 0
                row_idx = sum(isnan(H_DD), 2) == (self.sig_len - pg_num);
                H_DD(row_idx, :) = [];
            end
        end

        %{
        build the ideal pulse DD channel (callable after modulate)
        @taps_num:  the number of paths
        @his:       the channel gains
        @lis:       the channel delays
        @kis:       the channel dopplers
        %}
        function H_DD = buildIdealChannel(self, p, his, lis, kis)
            % input check
            if self.pulse_type ~= self.PULSE_IDEAL
                error("Cannot build the ideal pulse DD channel while not using ideal pulse.");
            end
            hw = zeros(self.nTimeslotNum, self.nSubcarNum);
            H_DD = zeros(self.sig_len, self.sig_len);
            for l = 1:self.nSubcarNum
                for k = 1:self.nTimeslotNum
                    for tap_id = 1:p
                        hi = his(tap_id);
                        li = lis(tap_id);
                        ki = kis(tap_id);
                        hw(k, l)= hw(k, l) + 1/self.sig_len*hi*exp(-2j*pi*li*ki/self.sig_len)*sum(exp(2j*pi*(l-li)*(0:self.nSubcarNum-1)/self.nSubcarNum))*sum(exp(-2j*pi*(k-ki)*(0:self.nTimeslotNum-1)/self.nTimeslotNum));
                    end
                    H_DD = H_DD + hw(k, l)*kron(circshift(eye(self.nTimeslotNum), k), circshift(eye(self.nSubcarNum), l));
                end
            end
        end

        %{
        build the rectangular pulse DD channel (callable after modulate)
        @taps_num:  the number of paths
        @his:       the channel gains
        @lis:       the channel delays
        @kis:       the channel dopplers
        %}
        function H_DD = buildRectaChannel(self, p, his, lis, kis)
            % input check
            if self.pulse_type ~= self.PULSE_RECTA
                error("Cannot build the rectangular pulse DD channel while not using rectanular pulse.");
            end
            % build H_DD
            H_DD = zeros(self.sig_len, self.sig_len); % intialize the return channel
            dftmat = dftmtx(self.nTimeslotNum)/sqrt(self.nTimeslotNum); % DFT matrix
            idftmat = conj(dftmat); % IDFT matrix 
            piMat = eye(self.sig_len); % permutation matrix (from the delay) -> pi
            % accumulate all paths
            for tap_id = 1:p
                hi = his(tap_id);
                li = lis(tap_id);
                ki = kis(tap_id);
                % delay
                piMati = circshift(piMat, li); 
                % Doppler
                %deltaMat_diag = exp(2j*pi*ki/(self.sig_len)*[0:self.sig_len-1-li, -li:-1]);
                deltaMat_diag = exp(2j*pi*ki/(self.sig_len)*[0:self.sig_len-1]);
                deltaMati = diag(deltaMat_diag);
                % Pi, Qi, & Ti
                Pi = kron(dftmat, eye(self.nSubcarNum))*piMati; 
                Qi = deltaMati*kron(idftmat, eye(self.nSubcarNum));
                Ti = Pi*Qi;
                H_DD = H_DD + hi*Ti;
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
end
