classdef OTFS < handle
    % constants
    properties(Constant)
        % Pilot
        PILOT_NO = 0;                               % no pilot
        PILOT_SINGLE_SISO = 10;                     % a single pilot for SISO case
        PILOT_MULTIP_SISO = 20;                     % multiple pilots for SISO case
        PILOT_TYPES = [OTFS.PILOT_NO, OTFS.PILOT_SINGLE_SISO, OTFS.PILOT_MULTIP_SISO];
        % Pilot locations
        PILOT_LOC_CENTER = 10;                      % the pilot is put at the center of frame
        PILOT_LOC_DELAY_MOST_CENTER = 20;           % the pilot is put at the most delay area center
        PILOT_LOC_TYPES = [OTFS.PILOT_LOC_CENTER, OTFS.PILOT_LOC_DELAY_MOST_CENTER];
        % Guard
        GUARD_NO = 0;
        GUARD_ENTIRE = 10;                          % guards take all rest positions
        GUARD_FULL = 20;                            % guards take all rows along the Doppler(time) axis
        GUARD_REDUCED = 30;                         % guards take limited rows along the Doppler(time) axis
        GUARD_TYPES = [OTFS.GUARD_NO, OTFS.GUARD_ENTIRE, OTFS.GUARD_FULL, OTFS.GUARD_REDUCED];
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
        y_DD                                        % Rx value in the DD domain (vectorized)
        taps_num = 0                                % paths number              
        delay_taps                                  % delay index, a row vector
        doppler_taps                                % doppler index (integers or fractional numbers), a row vector
        chan_coef                                   % path gain, a row vector
        % invalid area in X_DD (these parameters will be set after `insertPilotsAndGuards`)
        X_DD_invalid_num = 0;
        X_DD_invalid_delay_beg = NaN;
        X_DD_invalid_delay_end = NaN;
        X_DD_invalid_doppl_beg = NaN;
        X_DD_invalid_doppl_end = NaN;
        % pilot
        pilot_type = OTFS.PILOT_NO;                 % pilot - type
        pilots = [];                                % pilots value
        pilots_num_delay = 0;                       % pilots number along the delay(Doppler) axis
        pilots_num_doppler = 0;                     % pilots number along the Doppler(time) axis
        % pilot location
        pilot_loc_type = OTFS.PILOT_LOC_CENTER;
        % guard
        guard_type = OTFS.GUARD_NO;                 % guard - type
        % channel estimation
        ce_delay_beg = NaN;
        ce_delay_end = NaN;
        ce_doppl_beg = NaN;
        ce_doppl_end = NaN;
        ce_delay_taps                               % estimated delay index, a row vector
        ce_doppler_taps                             % estimated doppler index (integers or fractional numbers), a row vector
        ce_chan_coef                                % estimated path gain, a row vector
    end
    
    methods
        %% General OTFS Methods
        % constructor
        % @nSubcarNum:              subcarrier number
        % @nTimeslotNum:            timeslot number
        % @freq_spacing:            frequency spacing (kHz), the default is 15kHz
        % @fc:                      single carrier frequency (GHz), the default is 3GHz
        % @pilot_type:              pilot type
        % @pilot_loc_type:          pilot location type
        % @guard_type:              guard type
        function self = OTFS(nSubcarNum, nTimeslotNum, varargin)
            % Inputs Name-Value Pair 
            inPar = inputParser;
            addParameter(inPar, 'freq_spacing', self.freq_spacing, @isnumeric);   % register "freq_spacing"
            addParameter(inPar, 'fc', self.fc, @isnumeric);                       % register "fc"
            addParameter(inPar, "pilot_type", self.pilot_type, @(x) isscalar(x)&&isnumeric(x)&&ismember(x, self.PILOT_TYPES));
            addParameter(inPar, "pilot_loc_type", self.pilot_loc_type, @(x) isscalar(x)&&isnumeric(x)&&ismember(x, self.PILOT_LOC_TYPES));
            addParameter(inPar, "guard_type", self.guard_type, @(x) isscalar(x)&&isnumeric(x)&&ismember(x, self.GUARD_TYPES));
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
            % pillot
            self.pilot_type = inPar.Results.pilot_type;
            % guard
            self.guard_type = inPar.Results.guard_type;
            if self.pilot_type == self.PILOT_NO && self.guard_type ~= self.GUARD_NO
                error("Cannot allocate guards while not using pilots.");
            end
        end

        % modulate
        % @symbols: a vector of symbols to send or a matrix of [Doppler, delay] or [nTimeslotNum ,nSubcarNum]
        function modulate(self, symbols)
            % input check
            % input check - self.X_DD_invalid_num
            if self.pilot_type ~= self.PILOT_NO && self.X_DD_invalid_num == 0
                error("The modulation requires to insert pilots and guards first.");
            end
            data_num = self.nTimeslotNum*self.nSubcarNum - self.X_DD_invalid_num;
            % input check - symbols
            if ~isvector(symbols) && ~ismatrix(symbols)
                error("The transmission symbol must be a vector or a matrix");
            elseif isvector(symbols)
                % input check - symbols - vector
                if length(symbols) ~= data_num
                    error("The transmission symbol number must be %d", data_num);
                end
            elseif ismatrix(symbols) 
                % input check - symbols - matrix
                [sym_row_num, sym_col_num] = size(symbols);
                if sym_row_num ~= self.nTimeslotNum || sym_col_num ~= self.nSubcarNum
                    error("The transmission symbol must be in the shape of (%d, %d)", self.nTimeslotNum, self.nSubcarNum);
                end
                if data_num ~= self.nTimeslotNum*self.nSubcarNum
                    if sum(abs(symbols(self.X_DD_invalid_doppl_beg:self.X_DD_invalid_doppl_end, self.X_DD_invalid_delay_beg:self.X_DD_invalid_delay_end)) >= eps) ~= 0
                        error("The transmission symbols[%d:%d, %d:%d] must be zero", self.X_DD_invalid_doppl_beg, self.X_DD_invalid_doppl_end, self.X_DD_invalid_delay_beg, self.X_DD_invalid_delay_end);
                    end
                end
            end
            
            % modulate
            % reshape(rowwise) to [Doppler, delay] or [nTimeslotNum ,nSubcarNum]
            if isvector(symbols)
                symbols = symbols(:);
                if data_num == self.nSubcarNum*self.nTimeslotNum
                    self.X_DD = transpose(reshape(symbols, self.nSubcarNum, self.nTimeslotNum));
                else
                    symbols_id = 1;
                    for doppl_id = 1:self.nTimeslotNum
                        for delay_id = 1:self.nSubcarNum
                            if doppl_id<self.X_DD_invalid_doppl_beg || doppl_id>self.X_DD_invalid_doppl_end || delay_id<self.X_DD_invalid_delay_beg || delay_id>self.X_DD_invalid_delay_end
                                self.X_DD(doppl_id, delay_id) = symbols(symbols_id);
                                symbols_id = symbols_id + 1;
                            end
                        end
                    end
                    assert(symbols_id - 1 == data_num);
                end
            else
                if data_num == self.nSubcarNum*self.nTimeslotNum
                    self.X_DD = symbols;
                else
                    for doppl_id = 1:self.nTimeslotNum
                        for delay_id = 1:self.nSubcarNum
                            if ~ismember(doppl_id, self.X_DD_invalid_doppl_beg:self.X_DD_invalid_doppl_end) && ~ismember(delay_id, self.X_DD_invalid_delay_beg:self.X_DD_invalid_delay_end)
                                self.X_DD(doppl_id, delay_id) = symbols(doppl_id, delay_id);
                            end
                        end
                    end
                end
            end
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
        
        %% Channel estimation
        % insert pilots and guards
        % @pilots:                  a vector of your pilots (if given `pilots_pow` won't be used)
        % @pilots_pow:              pilot power to generate random pilots
        % @pilots_num_delay:        pilots number along the delay(Doppler) axis
        % @pilots_num_doppler:      pilots number along the Doppler(time) axis
        % @guard_delay_num_neg:     guard number negatively along the delay(frequency) axis
        % @guard_delay_num_pos:     guard number positively along the delay(frequency) axis
        % @guard_doppler_num_neg:   guard number negatively along the Doppler(time) axis
        % @guard_doppler_num_pos:   guard number positively along the Doppler(time) axis
        function insertPilotsAndGuards(self, pilots_num_delay, pilots_num_doppler, varargin)
            % optional inputs - register
            inPar = inputParser;
            addParameter(inPar, 'pilots', self.pilots, @(x) isvector(x)&&isnumeric(x));
            addParameter(inPar, 'pilots_pow', NaN, @(x) isscalar(x)&&isnumeric(x));
            addParameter(inPar, 'guard_delay_num_neg', 0, @(x) isscalar(x)&&isnumeric(x));
            addParameter(inPar, 'guard_delay_num_pos', 0, @(x) isscalar(x)&&isnumeric(x));
            addParameter(inPar, 'guard_doppler_num_neg', 0, @(x) isscalar(x)&&isnumeric(x));
            addParameter(inPar, 'guard_doppler_num_pos', 0, @(x) isscalar(x)&&isnumeric(x));
            inPar.KeepUnmatched = true;     % Allow unmatched cases
            inPar.CaseSensitive = false;    % Allow capital or small characters
            parse(inPar, varargin{:});
            % optional inputs - assign
            self.pilots = inPar.Results.pilots;
            pilots_len = length(self.pilots);
            pilots_pow = inPar.Results.pilots_pow;
            guard_delay_num_neg = inPar.Results.guard_delay_num_neg;
            guard_delay_num_pos = inPar.Results.guard_delay_num_pos;
            guard_doppler_num_neg = inPar.Results.guard_doppler_num_neg;
            guard_doppler_num_pos = inPar.Results.guard_doppler_num_pos;
            % input check
            % input check - pilots
            if self.pilot_type == self.PILOT_NO
                % overwrite pilot settings
                pilots_num_delay = 0;
                pilots_num_doppler = 0;
                self.pilots = [];
            else
                if pilots_num_delay <= 0 || pilots_num_doppler <= 0
                    error("Pilot number on each axis must be positive");
                elseif ~isempty(self.pilots)
                    if pilots_len ~= pilots_num_delay*pilots_num_doppler
                        error("The manual pilot input do not have the required number.");
                    elseif pilots_len > self.nTimeslotNum*self.nSubcarNum
                        error("The manual pilot input overflows (over the OTFS frame size).");
                    end
                end
            end
            % input check - pilot power
            if pilots_len == 0 && isnan(pilots_pow)
                error("The pilots (linear) power is required while no manual pilot input.");
            end
            % input check - guard
            switch self.guard_type
                case self.GUARD_NO
                    % overwrite guard settings
                    guard_delay_num_neg = 0;
                    guard_delay_num_pos = 0;
                    guard_doppler_num_neg = 0;
                    guard_doppler_num_pos = 0;
                case self.GUARD_ENTIRE
                    % overwrite guard settings
                    guard_delay_num_neg = floor((self.nSubcarNum - pilots_num_delay)/2);
                    guard_delay_num_pos = self.nSubcarNum - pilots_num_delay - guard_delay_num_neg;
                    guard_doppler_num_neg = floor((self.nTimeslotNum - pilots_num_doppler)/2);
                    guard_doppler_num_pos = self.nTimeslotNum - pilots_num_doppler - guard_doppler_num_neg;
                case self.GUARD_FULL
                    % overwrite guard settings
                    guard_doppler_num_neg = floor((self.nTimeslotNum - pilots_num_doppler)/2);
                    guard_doppler_num_pos = self.nTimeslotNum - pilots_num_doppler - guard_doppler_num_neg;
                case self.GUARD_REDUCED
                    if guard_delay_num_neg <= 0 ||  guard_delay_num_pos <= 0 || guard_doppler_num_neg <= 0 || guard_doppler_num_pos <=0
                        error("Guard number along the delay/Doppler axis must be positive.");
                    end
                    if floor(guard_delay_num_neg) ~= guard_delay_num_neg ||  floor(guard_delay_num_pos) ~= guard_delay_num_pos ||  floor(guard_doppler_num_neg) ~= guard_doppler_num_neg || floor(guard_doppler_num_pos) ~= guard_doppler_num_pos
                        error("Guard number along the delay/Doppler axis must be integers.");
                    end
            end

            % initiate X_DD if empty
            if isempty(self.X_DD)
                self.X_DD = zeros(self.nTimeslotNum, self.nSubcarNum);
            end
            % initiate pilots if empty
            if isempty(self.pilots)
                pilots_len = pilots_num_delay*pilots_num_doppler;
                switch self.pilot_type
                    % SISO cases (only use ones)
                    case {self.PILOT_SINGLE_SISO, self.PILOT_MULTIP_SISO}
                        self.pilots = sqrt(pilots_pow/2)*(1+1j)*ones(pilots_len, 1);
                end
            end
            % allocate pilots
            if pilots_len == 0
                % no pilots no operation
                self.X_DD_invalid_num = 0;
            else
                % some pilots
                % calulate the data number for two axises
                data_delay_num = (self.nSubcarNum - pilots_num_delay - guard_delay_num_neg - guard_delay_num_pos);
                data_doppler_num = self.nTimeslotNum - pilots_num_doppler - guard_doppler_num_neg - guard_doppler_num_pos;
                % calculate the pilot start point shift due to asymmetric guards
                pilot_shift_delay_pos = guard_delay_num_pos - guard_delay_num_neg; 
                pilot_shift_doppler_pos = guard_doppler_num_pos - guard_doppler_num_neg;
                % locate the 1st coordinates of the pilots (following the positive direction of axises)
                plots_loc_delay = 0;
                plots_loc_doppler = 0;
                switch self.pilot_loc_type
                    case self.PILOT_LOC_CENTER
                        plots_loc_delay = floor(data_delay_num/2) + guard_delay_num_neg + pilot_shift_delay_pos + 1;
                        plots_loc_doppler = floor(data_doppler_num/2) + guard_doppler_num_neg + pilot_shift_doppler_pos + 1;
                    case self.PILOT_LOC_DELAY_MOST_CENTER 
                        plots_loc_delay = floor(data_delay_num/2) + guard_delay_num_neg + pilot_shift_delay_pos + 1;
                        plots_loc_doppler = data_doppler_num + guard_doppler_num_neg + 1;
                end
                % allocate pilots
                self.X_DD(plots_loc_doppler:plots_loc_doppler+pilots_num_doppler-1, plots_loc_delay:plots_loc_delay+pilots_num_delay-1) = transpose(reshape(self.pilots, pilots_num_delay, pilots_num_doppler));
                % calculate the invalid area in X_DD
                self.X_DD_invalid_num = (pilots_num_delay+guard_delay_num_neg+guard_delay_num_pos)*(pilots_num_doppler+guard_doppler_num_neg+guard_doppler_num_pos);
                self.X_DD_invalid_delay_beg = plots_loc_delay - guard_delay_num_neg;
                self.X_DD_invalid_delay_end = plots_loc_delay + pilots_num_delay - 1 + guard_delay_num_pos;
                self.X_DD_invalid_doppl_beg = plots_loc_doppler - guard_doppler_num_neg;
                self.X_DD_invalid_doppl_end = plots_loc_doppler + pilots_num_doppler - 1 + guard_doppler_num_pos;
                assert(self.X_DD_invalid_num == (self.X_DD_invalid_delay_end - self.X_DD_invalid_delay_beg + 1)*(self.X_DD_invalid_doppl_end - self.X_DD_invalid_doppl_beg + 1));
            end
            % calulate channel estimate area
            if pilots_len > 0
                switch self.guard_type
                    case self.GUARD_ENTIRE
                        self.ce_delay_beg = 1;
                        self.ce_delay_end = self.nSubcarNum;
                        self.ce_doppl_beg = 1;
                        self.ce_doppl_end = self.nTimeslotNum;
                    case self.GUARD_FULL
                        self.ce_delay_beg = self.X_DD_invalid_delay_beg + guard_delay_num_neg;
                        self.ce_delay_end = self.X_DD_invalid_delay_end;
                        self.ce_doppl_beg = 1;
                        self.ce_doppl_end = self.nTimeslotNum;
                    case self.GUARD_REDUCED
                        self.ce_delay_beg = self.X_DD_invalid_delay_beg + guard_delay_num_neg;
                        self.ce_delay_end = self.X_DD_invalid_delay_end;
                        self.ce_doppl_beg = self.X_DD_invalid_doppl_beg + floor(guard_doppler_num_neg/2);
                        self.ce_doppl_end = self.X_DD_invalid_doppl_end - floor(guard_doppler_num_pos/2);
                end
            end
        end
        
        % estimated the channel
        % @threshold: the 
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
            % input check - pilot
            if self.pilot_type == self.PILOT_NO
                error("The channel estimation cannot be done because no pilot is in use.");
            end
            % input check - insert pilots & guards
            if isnan(self.ce_delay_beg) || isnan(self.ce_delay_end) || isnan(self.ce_doppl_beg) || isnan(self.ce_doppl_end)
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
            gains = [];
            delays = [];
            dopplers = [];
            switch self.pilot_type
                case self.PILOT_SINGLE_SISO
                    for delay_id = self.ce_delay_beg:self.ce_delay_end
                        for doppl_id = self.ce_doppl_beg:self.ce_doppl_end
                            tap_value = self.Y_DD(doppl_id, delay_id);
                            if abs(tap_value) > threshold
                                gains(end+1) = tap_value/self.pilots/exp(1j*2*pi*( - 1)*(k_id - pilot_loc_N)/M/N)
                            sfsfs/ewew?>.43
                            end
                        end
                    end
                case self.PILOT_MULTIP_SISO
                    error("Multiple pilots CE is not supported.");
            end
            
            
        end
        
        %% OTFS Detectors
        % detect
        % @detect_type:             detect type
        % @csi_type:                detect CSI type
        % @No:                      estimated noise power (linear)
        % @constellation:           the constellation (a vector)
        % @chan_coef:               the estimated channel gains
        % @delay_taps:              the estimated channel delays
        % @doppler_taps:            the estimated channel 
        function symbols = detect(self, detect_type, csi_type, No, constellation, varargin)
            % optional inputs - register
            inPar = inputParser;
            addParameter(inPar,"chan_coef", [], @isnumeric);
            addParameter(inPar,"delay_taps", [], @isnumeric);
            addParameter(inPar,"doppler_taps", [], @isnumeric);
            inPar.KeepUnmatched = true;
            inPar.CaseSensitive = false;
            parse(inPar, varargin{:});
            % optional inputs - assign
            chan_coef_est = inPar.Results.chan_coef;
            delay_taps_est = inPar.Results.delay_taps;
            doppler_taps_est = inPar.Results.doppler_taps;
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
                if self.pilot_type == self.PILOT_NO
                    error("Need CSI from channel estimation but no pilot is assigned.");
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
                constellation = constellation(:);
                constellation = constellation.';
            end
            
            % retrieve the channel information
            switch csi_type
                case self.DETECT_CSI_PERFECT
                    chan_coef_est = self.chan_coef;
                    delay_taps_est = self.delay_taps;
                    doppler_taps_est = self.doppler_taps;
                case self.DETECT_CSI_CE
                    
            end
            % estimate the symbols
            switch detect_type
                case self.DETECT_MP_BASE
                    symbols = self.detectMPBase(No, constellation, chan_coef_est, delay_taps_est, doppler_taps_est);
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
        % get the channel state information
        % @sort_by_gain: sort axis
        % @sort_by_delay_doppler: sort axes
        % @sort_by_doppler_delay: sort axes
        % @descend: sort direction
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
        
        % sort arrays based on the 
        % @arrs:                arrays [arr_num, arr_data_num]
        % @first_order_arr:     the 1st array id we follow the order
        % @second_order_arr:    the 2nd array id we follow the order
        % @descend:             descend order
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
            arr_data_id = 1;
            while arr_data_id <= arr_data_num
                % find the beginning and the end of the segment
                seg_beg = find(arr1, arr1(arr_data_id), "first");
                seg_end = find(arr1, arr1(arr_data_id), "last");
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
                    segs = segs(arr_id, :);
                    segs(arr_id, :) = segs(sort_idx_dim2);
                end
                % arrs - sort
                arrs(:, seg_beg:seg_end) = segs;
            end
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
