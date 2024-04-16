classdef OTFSDetector < handle
    properties(Constant)
        % Detect
        DETECT_MP_BASE = 10;                            % base OTFS MP detector proposed by P. Raviteja in 2018
        DETECT_TYPES = [OTFSDetector.DETECT_MP_BASE];
    end
    properties
        nSubcarNum {mustBeInteger}          % subcarrier number
        nTimeslotNum {mustBeInteger}        % timeslot number
        constellation                       % constellation
        constellation_len = 0;
        detect_type;
        x_num
        y
        y_num
        HDD
        p                                   % path number
        his
        lis
        kis
        No = NaN;
        % detectors
        % detector - mp base
        mp_base_n_ite = 200;
        mp_base_delta_fra = 0.6;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % General Methods
    methods
        %{
        constructor
        @nSubcarNum:                        subcarrier number
        @nTimeslotNum:                      timeslot number
        @constellation:                     the constellation (a vector)
        %}
        function self = OTFSDetector(nSubcarNum, nTimeslotNum, constellation)
            % nSubcarNum & nTimeslotNum
            if nSubcarNum <= 0
                error("Subcarrier number must be positive.");
            elseif nSubcarNum ~= floor(nSubcarNum)
                error("Subcarrier number must be an integer.");
            end
            self.nSubcarNum = nSubcarNum;
            if nTimeslotNum <= 0
                error("Timeslot number must be positive.");
            elseif nTimeslotNum ~= floor(nTimeslotNum)
                error("Timeslot number must be an integer.");
            end
            self.nTimeslotNum = nTimeslotNum;
            % constellation
            if ~isvector(constellation)
                error("The constellation must be a vector.");
            else
                % to row vector
                self.constellation = constellation(:).';
                self.constellation_len = length(constellation);
            end
        end
        
        %{
        detect
        @y:             the received signal from the resource grid (after demapping)
        @csi_info1:     (1) a vector of path gains  (2) a matrix of HDD
        @csi_info2:     (1) the delay indices       (2) the noise power
        @csi_info3:     (1) the Doppler indices
        @csi_info4:     (1) the noise power
        %}
        function symbols = detect(self, y, csi_info1, varargin)
            % input
            % input - y
            if ~isvector(y)
                error("The received signal from the resource grid (after demapping) must be a vector.");
            end
            self.y = y;
            self.y_num = length(y);
            csi_info_extra_len = length(varargin);
            % input - csi
            if isvector(csi_info1)
                % csi - type 1
                self.p = length(csi_info1);
                self.his = csi_info1;
                if csi_info_extra_len < 2
                    error("Still need at least two more inputs (the delay indices and the Doppler indices).");
                else
                    self.lis = varargin{1};
                    self.kis = varargin{2};
                    if self.p ~= length(self.lis) || self.p ~= length(self.kis)
                        error("The CSI must have the same length.");
                    end
                end
                if csi_info_extra_len >= 3
                    self.No = varargin{3};
                end
            elseif ismatrix(csi_info1)
                % csi - type 2
                self.HDD = csi_info1;
                [y_num_HDD, self.x_num] = size(self.HDD);
                if y_num_HDD ~= self.y_num
                    error("The HDD row number does not equal to the the received signal number.");
                end
                if csi_info_extra_len >= 1
                    self.No = varargin{1};
                end
            else
                error("The 1st CSI parameter must be a vector of path gains or a matrix of HDD.");
            end
            % input - csi - No
            if ~isnan(self.No)
                if ~isscalar(self.No)
                    error("The noise power must be a scalar.");
                end
            end
            % detect
            % detect - mp base
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Getters & Setters
    methods
        %{
        set detector types - MP Base
        @n_ite:           the iteration number (200 by default)
        @delta_fra:       the percentage for taking the values in the current iteration
        %}
        function useMPBase(self, varargin)
            self.detect_type = self.DETECT_MP_BASE; 
            % MP Base settings
            inPar = inputParser;
            addParameter(inPar,"n_ite", self.mp_base_n_ite, @(x) isscalar(x)&&isnumeric(x));
            addParameter(inPar,"delta_fra", self.mp_base_delta_fra, @(x) isscalar(x)&&isnumeric(x));
            inPar.KeepUnmatched = true;
            inPar.CaseSensitive = false;
            parse(inPar, varargin{:});
            self.mp_base_n_ite = inPar.Results.n_ite;
            self.mp_base_delta_fra = inPar.Results.delta_fra;
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % detectors
    methods(Access=private)
        %{
        MP base (proposed by P. Raviteja in 2017) from Emanuele Viterbo Research Group
        @No:              the estimated noise power
        @chan_coef:       the channel coefficient
        @delay_taps:      the delay indices
        @Doppler_taps:    the Doppler indices
        %}
        function symbols = detectMPBase(self, No, chan_coef, delay_taps, Doppler_taps, varargin)
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
end