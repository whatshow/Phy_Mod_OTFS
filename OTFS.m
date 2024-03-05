classdef OTFS < handle
    % constants
    properties(Constant)
        % Pilot
        PILOT_NO = 0;                               % no pilot
        PILOT_SINGLE_SISO = 1;                      % a single pilot for SISO case
        PILOT_TYPES = [OTFS.PILOT_NO, OTFS.PILOT_SINGLE_SISO];
        % Detect
        DETECT_NO = 0;                              % no detection
        DETECT_MP_BASE = 1;                         % base OTFS MP detector proposed by P. Raviteja in 2018
        DETECT_TYPES = [OTFS.DETECT_NO, OTFS.DETECT_MP_BASE];
        % CSI
        DETECT_CSI_PERFECT = 1;                     % perfect CSI
        DETECT_CSI_CE = 2;                          % CSI from channel estimation (its type is based on pilot type)
        DETECT_CSI_TYPES = [OTFS.DETECT_CSI_PERFECT, OTFS.DETECT_CSI_CE];
    end
    properties
        nSubcarNum {mustBeInteger}                  % subcarrier number
        nTimeslotNum {mustBeInteger}                % timeslot number
        freq_spacing {mustBeInteger} = 15           % frequency spacing (kHz), the default is 15kHz
        fc {mustBeInteger} = 3                      % single carrier frequency (GHz), the default is 3GHz
        X_DD                                        % Tx value in the delay Doppler(DD) domain
        X_TF                                        % Tx value in the time-frequency(TF) domain
        s                                           % Tx value in the time domain (array)
        H                                           % channel in the time domain
        r                                           % Rx value in the time domain (array)
        Y_TF                                        % Rx value in the TF domain
        Y_DD                                        % Rx value in the DD domain
        taps_num = 0                                % paths number              
        delay_taps                                  % delay index, a row vector
        doppler_taps                                % doppler index (integers or fractional numbers), a row vector
        chan_coef                                   % path gain, a row vector
        pilot_type = OTFS.PILOT_NO;                 % pilot - type 
        detect_type = OTFS.DETECT_NO;               % detect - type
        detect_csi_type = OTFS.DETECT_CSI_PERFECT;  % detect - CSI type
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
            addParameter(inPar,"pilot_type", self.pilot_type, @(x) isscalar(x)&&isnumeric(x)&&ismember(x, OTFS.PILOT_TYPES));
            addParameter(inPar,"detect_type", self.detect_type, @(x) isscalar(x)&&isnumeric(x)&&ismember(x, OTFS.DETECT_TYPES));
            addParameter(inPar,"detect_csi_type", self.detect_csi_type, @(x) isscalar(x)&&isnumeric(x)&&ismember(x, OTFS.DETECT_CSI_TYPES));
            inPar.KeepUnmatched = true;                                          % Allow unmatched cases
            inPar.CaseSensitive = false;                                         % Allow capital or small characters
            
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
            % freq_spacing & fc
            parse(inPar, varargin{:}); 
            self.freq_spacing = inPar.Results.freq_spacing;
            self.fc = inPar.Results.fc;
            % detect, pillot & CSI
            self.pilot_type = inPar.Results.pilot_type;
            self.detect_type = inPar.Results.detect_type;
            self.detect_csi_type = inPar.Results.detect_csi_type;
            if self.detect_type ~= OTFS.DETECT_NO
                if self.detect_csi_type == OTFS.DETECT_CSI_CE && self.pilot_type == OTFS.PILOT_NO
                    error("Cannot detect symbols while not use pilots but need CSI from channel estiamtion.");
                end
            end
        end

        % modulate
        % @symbols: a vector of symbols to send or a matrix of [Doppler, delay] or [nTimeslotNum ,nSubcarNum]
        function modulate(self, symbols)
            % input check
            if isvector(symbols)
                if length(symbols) ~= self.nTimeslotNum*self.nSubcarNum
                    error("The transmission symbol number must be %d", self.nTimeslotNum*self.nSubcarNum);
                end
            elseif ismatrix(symbols) 
                [sym_row_num, sym_col_num] = size(symbols);
                if sym_row_num ~= self.nTimeslotNum || sym_col_num ~= self.nSubcarNum
                    error("The transmission symbol must be in the shape of (%d, %d)", self.nTimeslotNum, self.nSubcarNum);
                end
            else
                error("The transmission symbol must be a vector or a matrix");
            end
            
            % modulate
            % reshape(rowwise) to [Doppler, delay] or [nTimeslotNum ,nSubcarNum]
            if isvector(symbols)
                symbols = symbols(:);
                self.X_DD = transpose(reshape(symbols, self.nSubcarNum, self.nTimeslotNum));
            else
                self.X_DD = symbols;
            end
            %self.X_DD = reshape(symbols, self.nTimeslotNum, self.nSubcarNum);
            % ISFFT 
            X_FT = fft(ifft(self.X_DD).').'/sqrt(self.nSubcarNum/self.nTimeslotNum);
            % X_TF is [nSubcarNum, nTimeslotNum]
            self.X_TF = X_FT.'; 
            % Heisenberg transform
            s_mat = ifft(self.X_TF)*sqrt(self.nSubcarNum);
            % vectorize
            self.s = s_mat(:);
        end
        
        % demodulate
        function yDD = demodulate(self)
            r_mat = reshape(self.r, self.nSubcarNum, self.nTimeslotNum);
            % Wigner transform (Y_TF in [nSubcarNum, nTimeslotNum])
            self.Y_TF = fft(r_mat)/sqrt(self.nSubcarNum); 
            Y_FT = self.Y_TF.';
            % SFFT (Y_DD in [Doppler, delay] or [nTimeslotNum ,nSubcarNum])
            self.Y_DD = ifft(fft(Y_FT).').'/sqrt(self.nTimeslotNum/self.nSubcarNum); 
            % return the DD domain vector
            yDD = self.Y_DD.';
            yDD = yDD(:);
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
        % @noisePow: noise power (a scalar)
        function r = passChannel(self, noPow)
            % input check
            if ~isscalar(noPow)
                error("The noise power must be a scalar.");
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
            % add noise
            if noPow >= 0
                noise = sqrt(noPow/2)*(randn(size(s_chan)) + 1i*randn(size(s_chan)));
                self.r = s_chan + noise;
            else
                self.r = s_chan;
            end
            % remove CP
            self.r = self.r(cp_len+1:cp_len+(self.nTimeslotNum*self.nSubcarNum));
            % return
            r = self.r;
        end
        
        %% OTFS Detectors
        % detect
        % @No: estimated noise power (linear)
        % @constellation: the constellation (a vector)
        function symbols = detect(self, No, constellation, varargin)
            % input check
            if ~isscalar(No)
                error("The noise power(linear) must be a scalar.");
            end
            if ~isvector(constellation)
                error("The constellation must be a vector.");
            else
                % to row vector
                constellation = constellation(:);
                constellation = constellation.';
            end
            % input check - optional
            inPar = inputParser;
            addParameter(inPar,"H_DD", [], @(x) isempty(x)||~isvector(x)&&ismatrix(x)&&isnumeric(x));
            addParameter(inPar,"chan_coef", [], @(x) isempty(x)||isscalar(x)&&isnumeric(x));
            addParameter(inPar,"delay_taps", [], @(x) isempty(x)||isscalar(x)&&isnumeric(x));
            addParameter(inPar,"doppler_taps", [], @(x) isempty(x)||isscalar(x)&&isnumeric(x));
            inPar.KeepUnmatched = true;
            inPar.CaseSensitive = false;
            parse(inPar, varargin{:});
            % take inputs
            H_DD = inPar.Results.H_DD;
            chan_coef = inPar.Results.chan_coef;
            delay_taps = inPar.Results.delay_taps;
            doppler_taps = inPar.Results.doppler_taps;
            
            % retrieve the channel information
            switch self.detect_csi_type
                case OTFS.DETECT_CSI_PERFECT
                    chan_coef = self.chan_coef;
                    delay_taps = self.delay_taps;
                    doppler_taps = self.doppler_taps;
                case OTFS.DETECT_CSI_CE
                    
            end
            % estimate the symbols
            switch self.detect_type
                case OTFS.DETECT_NO
                    symbols = [];
                case OTFS.DETECT_MP_BASE
                    symbols = self.detectMPBase(No, constellation, chan_coef, delay_taps, doppler_taps);
            end
        end
        
        % MP base (proposed by P. Raviteja in 2017) from Emanuele Viterbo Research Group
        % @No:              the estimated noise power
        % @constellation:   the constellation (a vector)
        % @chan_coef:       the channel coefficient
        % @delay_taps:      the delay indices
        % @Doppler_taps:    the Doppler indices
        % @n_ite:           the iteration number (200 by default)
        % @delta_fra:       the percentage for taking the values in the current iteration
        function symbols = detectMPBase(self, No, constellation, chan_coef, delay_taps, Doppler_taps, varargin)
            % input check
            constellation_len = length(constellation);
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
            % set initial values
            yv = reshape(self.Y_DD, self.nTimeslotNum*self.nSubcarNum, 1);
            mean_int = zeros(self.nTimeslotNum*self.nSubcarNum,taps);
            var_int = zeros(self.nTimeslotNum*self.nSubcarNum,taps);
            p_map = ones(self.nTimeslotNum*self.nSubcarNum,taps, constellation_len)*(1/constellation_len);
            % detect
            conv_rate_prev = -0.1;
            for ite=1:n_ite
                % Update mean and var
                for ele1=1:1:self.nSubcarNum
                    for ele2=1:1:self.nTimeslotNum
                        mean_int_hat = zeros(taps,1);
                        var_int_hat = zeros(taps,1);
                        for tap_no=1:taps
                            m = ele1-1-delay_taps(tap_no)+1;
                            add_term = exp(1i*2*(pi/self.nSubcarNum)*(m-1)*(Doppler_taps(tap_no)/self.nTimeslotNum));
                            add_term1 = 1;
                            if ele1-1<delay_taps(tap_no)
                                n = mod(ele2-1-Doppler_taps(tap_no),self.nTimeslotNum) + 1;
                                add_term1 = exp(-1i*2*pi*((n-1)/self.nTimeslotNum));
                            end
                            new_chan = add_term * (add_term1) * chan_coef(tap_no);

                            for i2=1:1:constellation_len
                                mean_int_hat(tap_no) = mean_int_hat(tap_no) + p_map(self.nTimeslotNum*(ele1-1)+ele2,tap_no,i2) * constellation(i2);
                                var_int_hat(tap_no) = var_int_hat(tap_no) + p_map(self.nTimeslotNum*(ele1-1)+ele2,tap_no,i2) * abs(constellation(i2))^2;
                            end
                            mean_int_hat(tap_no) = mean_int_hat(tap_no) * new_chan;
                            var_int_hat(tap_no) = var_int_hat(tap_no) * abs(new_chan)^2;
                            var_int_hat(tap_no) = var_int_hat(tap_no) - abs(mean_int_hat(tap_no))^2;
                        end

                        mean_int_sum = sum(mean_int_hat);
                        var_int_sum = sum(var_int_hat)+(No);

                        for tap_no=1:taps
                            mean_int(self.nTimeslotNum*(ele1-1)+ele2,tap_no) = mean_int_sum - mean_int_hat(tap_no);
                            var_int(self.nTimeslotNum*(ele1-1)+ele2,tap_no) = var_int_sum - var_int_hat(tap_no);
                        end

                    end
                end
                %% Update probabilities
                sum_prob_comp = zeros(self.nTimeslotNum*self.nSubcarNum, constellation_len);
                dum_eff_ele1 = zeros(taps,1);
                dum_eff_ele2 = zeros(taps,1);
                for ele1=1:1:self.nSubcarNum
                    for ele2=1:1:self.nTimeslotNum
                        dum_sum_prob = zeros(constellation_len,1);
                        log_te_var = zeros(taps,constellation_len);
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
                            new_chan = add_term * add_term1 * chan_coef(tap_no);

                            dum_eff_ele1(tap_no) = eff_ele1;
                            dum_eff_ele2(tap_no) = eff_ele2;
                            for i2=1:1:constellation_len
                                dum_sum_prob(i2) = abs(yv(self.nTimeslotNum*(eff_ele1-1)+eff_ele2)- mean_int(self.nTimeslotNum*(eff_ele1-1)+eff_ele2,tap_no) - new_chan * constellation(i2))^2;
                                dum_sum_prob(i2)= -(dum_sum_prob(i2)/var_int(self.nTimeslotNum*(eff_ele1-1)+eff_ele2,tap_no));
                            end
                            dum_sum = dum_sum_prob - max(dum_sum_prob);
                            dum1 = sum(exp(dum_sum));
                            log_te_var(tap_no,:) = dum_sum - log(dum1);
                        end
                        for i2=1:1:constellation_len
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
                            p_map(self.nTimeslotNum*(eff_ele1-1)+eff_ele2,tap_no,:) = (dum_sum/dum1)*delta_fra + (1-delta_fra)*reshape(p_map(self.nTimeslotNum*(eff_ele1-1)+eff_ele2,tap_no,:),1,constellation_len);
                        end

                    end
                end
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
                    X_DD_est(ele2,ele1) = constellation(pos);
                end
            end
            % extract symbols
            symbols = X_DD_est.';
            symbols = symbols(:);
        end
        
        %% support function
        % Get the channel matrix in Delay Doppler Domain (using the rectangular waveform)
        function H_DD = getChannel(self)
            % intialize the return channel
            H_DD = zeros(self.nTimeslotNum*self.nSubcarNum, self.nTimeslotNum*self.nSubcarNum);
            % DFT & IDFT matrix
            dftmat = dftmtx(self.nTimeslotNum)/sqrt(self.nTimeslotNum);
            idftmat = conj(dftmat);
            % permutation matrix (from the delay) -> pi
            piMat = eye(self.nTimeslotNum*self.nSubcarNum);
            % accumulate all paths
            for tap_id = 1:self.taps_num
                hi = self.chan_coef(tap_id);
                li = self.delay_taps(tap_id);
                ki = self.doppler_taps(tap_id);
                % delay
                piMati = circshift(piMat, li);
                % Doppler
                deltaMat_diag = exp(1j*2*pi*ki/(self.nTimeslotNum*self.nSubcarNum)*(0:1:self.nTimeslotNum*self.nSubcarNum-1));
                deltaMati = diag(deltaMat_diag);
                % Pi
                Pi = kron(dftmat, eye(self.nSubcarNum))*piMati;
                % Qi
                Qi = deltaMati*kron(idftmat, eye(self.nSubcarNum));
                % Ti
                Ti = Pi*Qi;
                % add this path
                H_DD = H_DD + hi*Ti;
            end
        end
        % get the channel delays
        function delays = getChannelDelays(self)
            delays = self.delay_taps;
        end
        % get the channel dopplers
        function dopplers = getChannelDopplers(self)
            dopplers = self.doppler_taps;
        end
        % get the channel gains
        function gains = getChannelGains(self)
            gains = self.chan_coef;
        end
        % get the signal in the Delay Time domain [delay, time]
        function X_DT = getX2DT(self)
            X_DT = ifft(self.X_DD).';
        end
        % get the signal in the TF domain
        function X_TF = getX2TF(self)
            X_TF = self.X_TF;
        end
        % get the signal in the time domain
        function s = getX2T(self, varargin)
            % Inputs Name-Value Pair 
            inPar = inputParser;
            addParameter(inPar,'fft_size', 0, @isnumeric);
            inPar.KeepUnmatched = true;                                          % Allow unmatched cases
            inPar.CaseSensitive = false;                                         % Allow capital or small characters
            % freq_spacing & fc
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
        % get the received signal in delay Doppler domain
        function Y_DD = getYDD(self)
            Y_DD = self.Y_DD;
        end
    end
end

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
