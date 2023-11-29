classdef OTFS < handle
    properties
        nSubcarNum {mustBeInteger}                  % subcarrier number
        nTimeslotNum {mustBeInteger}                % timeslot number
        freq_spacing {mustBeInteger} = 15           % frequency spacing (kHz), the default is 15kHz
        fc {mustBeInteger} = 3                      % single carrier frequency (GHz), the default is 3GHz
        X_DD                                        % Tx value in the delay Doppler(DD) domain
        X_TF                                        % Tx value in the time-frequency(TF) domain
        s                                           % Tx value in the time domain (array)
        H                                           % channel in the time domain
        H_DD                                        % Effective channel in the DD domain
        r                                           % Rx value in the time domain (array)
        Y_TF                                        % Rx value in the TF domain
        Y_DD                                        % Rx value in the DD domain
        taps_num                                    % paths number              
        delay_taps                                  % delay index
        doppler_taps                                % doppler index (it can be integers or fractional numbers)
        chan_coef                                   % path gain
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
        % @symbols: a vector of symbols to send
        function modulate(self, symbols)
            % input check
            if length(symbols) ~= self.nTimeslotNum*self.nSubcarNum
                error("The transmission symbol number must be %d", self.nTimeslotNum*self.nSubcarNum);
            end
            
            % modulate
            % reshape(rowwise) to [Doppler, delay] or [nTimeslotNum ,nSubcarNum]
            symbols = symbols(:);
            self.X_DD = transpose(reshape(symbols, self.nSubcarNum, self.nTimeslotNum)); 
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
        % @p: the path number
        % @lmax: the maxmimal delay index
        % @kmax: the maximal Doppler index
        % @is_fractional_doppler: whether we use the fractional Doppler (the default is false)
        function H_DD = setChannel(self, p, lmax, kmax, varargin)
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
            elseif p > (lmax + 1)*(2*kmax+1)
                error("The path number must be less than lmax*(2*kmax+1) = %d", (lmax + 1)*(2*kmax+1));
            end
            % is fractional doppler
            parse(inPar, varargin{:});
            is_fractional_doppler = inPar.Results.is_fractional_doppler;
            
            % generate paths
            lmin= 0;
            kmin = -kmax;
            taps_max = (kmax - kmin + 1)*(lmax - lmin + 1);
            % create delay options [lmin, lmin, lmin, lmin+1, lmin+1, lmin+1 ...]
            delay_taps_all = kron(lmin:lmax, ones(1, kmax - kmin + 1)); 
            % create Doppler options [kmin, kmin+1, kmin+2 ... kmax, kmin ...]
            Doppler_taps_all = repmat(kmin:kmax, 1, lmax - lmin + 1);
            % We select P paths from all possible paths; that is, we do the randperm(taps_max) and we choose the first P items
            taps_idx_chaotic = randperm(taps_max);
            taps_selected_idx = taps_idx_chaotic(1:p);
            % set delay & Doppler
            self.delay_taps = delay_taps_all(taps_selected_idx);
            self.doppler_taps = Doppler_taps_all(taps_selected_idx);
            % set the path gain
            self.chan_coef = sqrt(1/p)*(sqrt(1/2) * (randn(1, p)+1i*randn(1, p)));
            % record path number
            self.taps_num = p;
            
            % return the channel
            H_DD = self.getChannel();
        end
        
        % pass the channel
        % @noisePow: noise power (a scalar)
        function passChannel(self, noPow)
            % input check
            if ~isscalar(noPow)
                error("The noise power must be a scalar.");
            end
            % add CP
            cp_len = max(self.delay_taps);
            s_cp = [self.s(self.nTimeslotNum*self.nSubcarNum - cp_len + 1 : self.nTimeslotNum*self.nSubcarNum);self.s];
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
        end
        
        %% support function
        % Get the channel matrix (using the rectangular waveform)
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
            % record
            self.H_DD = H_DD;
        end
        
    end

end