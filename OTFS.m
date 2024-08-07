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
    end
    properties
        % system
        fc = 4;                                     % Carrier frequency (GHz)
        fq_sp = 15;                                 % subcarrier spacing (kHz)
        res_k = 0;                                  % resolution - Doppler (kHz)
        res_l = 0;                                  % resolution - delay (ms)
        % RG infomation
        nSubcarNum {mustBeInteger}                  % subcarrier number
        nTimeslotNum {mustBeInteger}                % timeslot number
        sig_len = 0;                                % the total signal length
        rg
        % OTFS process signal
        X_DD
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
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % General Methods
    methods
        %{
        constructor
        @fc:        Carrier frequency (GHz)
        @fq_sp:     subcarrier spacing (kHz)
        %}
        function self = OTFS(varargin)
            % optional inputs - register
            inPar = inputParser;
            addParameter(inPar,"fc", self.fc, @(x) isscalar(x)&isnumeric(x)); 
            addParameter(inPar,"fq_sp", self.fq_sp, @(x) isscalar(x)&isnumeric(x)); 
            inPar.KeepUnmatched = true;
            inPar.CaseSensitive = false;
            parse(inPar, varargin{:});
            self.fc = inPar.Results.fc;
            self.fq_sp = inPar.Results.fq_sp;
        end

        %{
        pulse settings
        %}
        function setPulse2Ideal(self)
            self.pulse_type = self.PULSE_IDEAL;
        end
        function setPulse2Recta(self)
            self.pulse_type = self.PULSE_RECTA;
        end

        %{
        modulate (use fast method by default)
        @in1:       an OTFS resource grid or a 2D matrix [(batch_size), Doppler, delay]
        @isFast:    DD domain -> TD domain (no X_TF) 
        %}
        function modulate(self, in1, varargin)
            % optional inputs - register
            inPar = inputParser;
            addParameter(inPar,"isFast", true, @(x) isscalar(x)&islogical(x)); 
            inPar.KeepUnmatched = true;
            inPar.CaseSensitive = false;
            parse(inPar, varargin{:});
            isFast = inPar.Results.isFast;
            % input check
            if isa(in1, 'OTFSResGrid')
                % load RG
                [self.nSubcarNum, self.nTimeslotNum] = in1.getContentSize();
                self.calcRes();
                self.sig_len = self.nSubcarNum*self.nTimeslotNum;
                if in1.isPulseIdeal()
                    self.pulse_type = self.PULSE_IDEAL;
                elseif in1.isPulseRecta()
                    self.pulse_type = self.PULSE_RECTA;
                else
                    error("Pulse shaping is not given in the resource grid.");
                end
                self.rg = in1.clone();
                self.X_DD = in1.getContent();
            elseif ismatrix(in1)
                [self.nTimeslotNum, self.nSubcarNum] = size(in1);
                self.sig_len = self.nSubcarNum*self.nTimeslotNum;
                self.X_DD = in1;
            else
                error("The input must be an OTFS resource grid or a 2D matrix [(batch_size), Doppler, delay].");
            end
            
            % modulate
            if self.pulse_type == self.PULSE_RECTA
                if isFast
                    s_mat = ifft(self.X_DD).'*sqrt(self.nTimeslotNum);
                    self.s = s_mat(:);
                else
                    X_FT = fft(ifft(self.X_DD).').'/sqrt(self.nSubcarNum/self.nTimeslotNum); % ISFFT 
                    self.X_TF = X_FT.'; % X_TF is [nSubcarNum, nTimeslotNum]
                    s_mat = ifft(self.X_TF)*sqrt(self.nSubcarNum); % Heisenberg transform
                    self.s = s_mat(:);
                end
            end
        end

        %{
        set channel (in1, in2, in3)
        set a fixed chanel (at least two paths, if you want to add one fixed path, call `setChannelExtra`)
        @in1->his:        the path gains
        @in2->lis:        the delays
        @in3->kis:        the doppler shifts
        set a random channel (overwritten the channel setting; use Rayleigh fading if not select channel model)
        @in1->p:          the path number
        @in2->lmax:       the maxmimal delay index
        @in3->kmax:       the maximal Doppler index (can be fractional)
        @force_frac:      use fractional Doppler (force)
        @isAWGN:          use awgn
        @isRician:        use Rician fading
        %}
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
                    if force_frac
                        kmax_frac = 0.5;
                    end
                end
                % input check
                if p > (lmax + 1)*(2*kmax+1)
                    error("The path number must be less than lmax*(2*kmax+1) = %d", (lmax + 1)*(2*kmax+1));
                end
                lmin= 1;
                kmin = -kmax;
                taps_max = (kmax - kmin + 1)*(lmax - lmin + 1);
                % create delay options [lmin, lmin, lmin, lmin+1, lmin+1, lmin+1 ...]
                l_combs = kron(lmin:lmax, ones(1, kmax - kmin + 1));
                % create Doppler options [kmin, kmin+1, kmin+2 ... kmax, kmin ...]
                k_combs = repmat(kmin:kmax, 1, lmax - lmin + 1);
                % select P paths from all possible paths
                taps_selected_idx = self.shufSelectTopNIdx(taps_max, p);
                % CSI - delay
                self.delay_taps = l_combs(taps_selected_idx);
                self.delay_taps(find(self.delay_taps == min(self.delay_taps), 1)) = 0;
                % CSI - doppler
                self.doppler_taps = k_combs(taps_selected_idx);
                % add fractional Doppler
                if force_frac
                    doppler_taps_k_max_pos_idx = self.doppler_taps == kmax;
                    doppler_taps_k_max_neg_idx = self.doppler_taps == -kmax;
                    doppler_taps_k_other_idx = abs(self.doppler_taps)~= kmax;
                    frac_range_max_pos = rand(1, p)*(kmax_frac+0.5)-0.5;
                    frac_range_max_neg = rand(1, p)*(-kmax_frac-0.5)+0.5;
                    frac_range_others = rand(1, p) - 0.5;
                    frac_range_all = frac_range_max_pos.*doppler_taps_k_max_pos_idx + frac_range_max_neg.*doppler_taps_k_max_neg_idx + frac_range_others.*doppler_taps_k_other_idx;
                    self.doppler_taps = self.doppler_taps + frac_range_all;
                end
                % CSI - others
                self.taps_num = p;
                if isAWGN
                    self.chan_coef = sqrt(1/2/p)*(ones(1, p)+1i*ones(1, p));
                elseif isRician
                    self.chan_coef = [sqrt(1/2/p), sqrt(1/2/p)*(randn(1, p-1)+1i*randn(1, p-1))];
                else
                    self.chan_coef = sqrt(1/2/p)*(randn(1, p)+1i*randn(1, p)); % use Rayleigh fading by default
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
                s_chan = H_DD*reshape(self.X_DD.', self.sig_len, 1);
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
        function ou1 = demodulate(self, varargin)
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
                    % Wigner transform (Y_TF in [nSubcarNum, nTimeslotNum])
                    self.Y_TF = fft(r_mat)/sqrt(self.nSubcarNum);
                    Y_FT = self.Y_TF.';
                    % SFFT (Y_DD in [Doppler, delay] or [nTimeslotNum ,nSubcarNum])
                    self.Y_DD = ifft(fft(Y_FT).').'/sqrt(self.nTimeslotNum/self.nSubcarNum);
                end
            end
            if isempty(self.rg)
                ou1 = self.Y_DD;
            else
                self.rg.setContent(self.Y_DD);
                ou1 = self.rg;
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Getters & Setters
    methods
        %{
        build Phi - the channel estimation matrix
        @X:     the Tx matrix in DD domain ([batch_size], doppler, delay)
        @lmax:  the maximal delay
        @kmax:  the maximal Doppler
        %}
        function Phi = buildPhi(self, X, lmax, kmax)
            pmax = (lmax+1)*(2*kmax+1);                 % the number of all possible paths
            lis = kron(0:lmax, ones(1, 2*kmax + 1));    % the delays on all possible paths
            kis = repmat(-kmax:kmax, 1, lmax+1);        % the dopplers on all possible paths
            Phi = zeros(self.sig_len, pmax);            % the return matrix
            for yk = 1:self.nTimeslotNum
                for yl = 1:self.nSubcarNum
                    Phi_ri = (yk - 1)*self.nSubcarNum + yl;      % row id in Phi
                    for p_id = 1:pmax
                        % path delay and doppler
                        li = lis(p_id);
                        ki = kis(p_id);
                        % x(k, l)_
                        xl = yl - li;
                        if yl-1 < li
                            xl = xl + self.nSubcarNum;
                        end
                        xk = mod(yk - 1 - ki, self.nTimeslotNum) + 1;
                        % exponential part (pss_beta)
                        if self.pulse_type == self.PULSE_IDEAL
                            pss_beta = exp(-2j*pi*li*ki/self.nSubcarNum/self.nTimeslotNum);
                        elseif self.pulse_type == self.PULSE_RECTA
                            pss_beta = exp(2j*pi*(yl - li - 1)*ki/self.nSubcarNum/self.nTimeslotNum); % here, you must use `yl-li-1` instead of `xl-1` or there will be an error
                            if yl-1 < li
                                pss_beta = pss_beta*exp(-2j*pi*(xk-1)/self.nTimeslotNum);
                            end
                        end
                        % assign value
                        Phi(Phi_ri, p_id) = X(xk, xl)*pss_beta;
                    end
                end
            end  
        end

        %{
        build Hdd
        @his: the path gains
        @lmax:  the maximal delay
        @kmax:  the maximal Doppler
        @thres(opt): the threshold of a path (default 0)
        %}
        function Hdd = buildHdd(self, his, lmax, kmax, varargin)
            % register optional inputs 
            inPar = inputParser;
            addParameter(inPar,"thres", 0, @(x) isscalar(x)&isnumeric(x));
            inPar.KeepUnmatched = true;
            inPar.CaseSensitive = false;
            parse(inPar, varargin{:});
            thres = inPar.Results.thres;
            % other inputs
            pmax = (lmax+1)*(2*kmax+1);                 % the number of all possible paths
            lis = kron(0:lmax, ones(1, 2*kmax + 1));    % the delays on all possible paths
            kis = repmat(-kmax:kmax, 1, lmax+1);        % the dopplers on all possible paths
            % filter the path gains
            his(abs(his) < abs(thres)) = 0;
            % build the channel in DD domain
            switch self.pulse_type
                case self.PULSE_IDEAL
                    Hdd = self.buildIdealChannel(pmax, his, lis, kis);
                case self.PULSE_RECTA
                    Hdd = self.buildRectaChannel(pmax, his, lis, kis);
                otherwise
                    error("The pulse type is not recognised.");
            end
        end

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
                    if p == 0
                        error("The input CSI is empty.");
                    elseif p ~= length(lis) && p ~= length(kis)
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
            if data_only && ~isempty(self.rg)
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
    % private methods
    methods(Access=private)
        %{
        calculate the resolution
        %}
        function calcRes(self)
            self.res_k = self.fq_sp/self.nTimeslotNum;
            self.res_l = 1/self.fq_sp/self.nTimeslotNum;
        end

        %{
        shuffle and select top n elements' indices 
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
                s_cp_t = (0:self.sig_len-1)';
            elseif self.cp_type == self.CP_ONE_FRAM
                s_cp = [self.s(end - self.cp_len + 1 : end);self.s];
                s_cp_t = (-self.cp_len:self.sig_len-1)';
            elseif self.cp_type == self.CP_ONE_FRAM_SUB
                s_mat = reshape(self.s, self.nSubcarNum, self.nTimeslotNum);
                s_mat = [s_mat(end - self.cp_len + 1 : end, :); s_mat];
                s_cp = s_mat(:);
                s_cp_t = repmat((-self.cp_len:self.nSubcarNum-1)', self.nTimeslotNum, 1);
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
                s_chan_mat = reshape(s_chan, [], self.nTimeslotNum);
                s_chan_rm_mat = s_chan_mat(self.cp_len+1:self.cp_len+self.nSubcarNum, :);
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
            if isempty(self.rg)
                error("The resource grid is not given.");
            end

            %TODO: remove zero padding area
            invalud_row = nan(1, self.sig_len);
            invalud_col = nan(self.sig_len, 1);
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
                col_idx = sum(isnan(H_DD)) == self.sig_len;
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
                        hw_add = 1/self.sig_len*hi*exp(-2j*pi*li*ki/self.sig_len)* ...
                                sum(exp(2j*pi*(l-li)*(0:self.nSubcarNum-1)/self.nSubcarNum))* ...
                                sum(exp(-2j*pi*(k-ki)*(0:self.nTimeslotNum-1)/self.nTimeslotNum));
                        hw(k, l)= hw(k, l) + hw_add;
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
            H_DD = zeros(self.sig_len, self.sig_len);                       % intialize the return channel
            dftmat = dftmtx(self.nTimeslotNum)/sqrt(self.nTimeslotNum);     % DFT matrix
            idftmat = conj(dftmat);                                         % IDFT matrix 
            piMat = eye(self.sig_len);                                      % permutation matrix (from the delay) -> pi
            % accumulate all paths
            for tap_id = 1:p
                hi = his(tap_id);
                li = lis(tap_id);
                ki = kis(tap_id);
                % delay
                piMati = circshift(piMat, li); 
                % Doppler
                deltaMat_diag = exp(2j*pi*ki/(self.sig_len)*self.buildTimeSequence(li));
                deltaMati = diag(deltaMat_diag);
                % Pi, Qi, & Ti
                Pi = kron(dftmat, eye(self.nSubcarNum))*piMati; 
                Qi = deltaMati*kron(idftmat, eye(self.nSubcarNum));
                Ti = Pi*Qi;
                H_DD = H_DD + hi*Ti;
            end
        end

        %{
        build the time sequence for the given delay
        %}
        function ts = buildTimeSequence(self, li)
            ts = [0:self.sig_len-1-li, -li:-1];
        end
    end
end
