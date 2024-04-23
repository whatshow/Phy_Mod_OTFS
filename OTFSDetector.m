classdef OTFSDetector < handle
    properties(Constant)
        % Detect
        DETECT_MP_BASE = 10;                            % base OTFS MP detector proposed by P. Raviteja in 2018
        DETECT_TYPES = [OTFSDetector.DETECT_MP_BASE];
    end
    properties
        constel                       % constel
        constel_len = 0;
        M {mustBeInteger}                   % subcarrier number
        N {mustBeInteger}                   % timeslot number
        detect_type;
        x_num
        y_rg
        y
        y_num
        % CSI
        HDD = NaN;
        p = 0;
        his = [];
        lis = [];
        kis = [];
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
        @constel:           the constellation (a vector)
        @
        %}
        function self = OTFSDetector(constel)
            % constel
            if ~isvector(constel)
                error("The constel must be a vector.");
            else
                % to row vector
                self.constel = constel(:).';
                self.constel_len = length(constel);
            end
        end
        
        %{
        detect
        @y:             the received signal from the resource grid (after demapping) or just a resource grid 
        @csi_info1:     (1) a vector of path gains  (2) a matrix of HDD
        @csi_info2:     (1) the delay indices       (2) the noise power
        @csi_info3:     (1) the Doppler indices
        @csi_info4:     (1) the noise power
        %}
        function syms = detect(self, y, csi_info1, varargin)
            % input
            % input - y
            if isa(y, "OTFSResGrid")
                self.y_rg = y;
                % nSubcarNum & nTimeslotNum
                [self.M, self.N] = y.getContentSize();
                if self.M <= 0
                    error("Subcarrier number must be positive.");
                elseif self.M ~= floor(self.M)
                    error("Subcarrier number must be an integer.");
                end
                if self.N <= 0
                    error("Timeslot number must be positive.");
                elseif self.N ~= floor(self.N)
                    error("Timeslot number must be an integer.");
                end
                self.y_num = self.M*self.N;
            elseif isvector(y)
                self.y = y;
                self.y_num = length(y);
            else
                error("y must be a vector from the resource grid (after demapping) or an OTFS resource grid.");
            end
            
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
            if self.detect_type == self.DETECT_MP_BASE
                syms = self.detectMPBase();
            end
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
        %}
        function syms = detectMPBase(self)
            % input check
            if isnan(self.No) || ~isscalar(self.No)
                error("The noise power must be a scalar.");
            end
            % set data
            conv_rate_prev = -0.1;
            Y_DD = self.y_rg.getContent();
            [pg_num, pg_l_beg, pg_l_end, pg_k_beg, pg_k_end] = self.y_rg.getAreaPG();   % PG  
            [ce_num, ce_l_beg, ce_l_end, ce_k_beg, ce_k_end] = self.y_rg.getAreaCE();   % CE
        
            % init all parameters
            % mean & variance (d, c), d=M*(l-1)+k, c=M*(l-li-1)+(k-ki)
            mu_dc = zeros(self.N*self.M, self.p);
            sigma2_dc = zeros(self.N*self.M, self.p);
            % probability (d, c), d=M*(l-1)+k, c=M*(l-1+li)+(k+ki)
            p_dc = ones(self.N*self.M, self.p, self.constel_len)*(1/self.constel_len);
            for ite=1:self.mp_base_n_ite
                % we update the mean and variance for each d in mu[d, c]
                for l=1:self.M
                    for k=1:self.N
                        d = self.N*(l-1)+k;
                        % jump if y[d] in CE
                        if ce_num>0 && l>=ce_l_beg && l<=ce_l_end && k>=ce_k_beg && k<=ce_k_end
                            continue;
                        end
                        mu_d = zeros(self.p,1);                         % the sum of mu[d, c] for a given d (must be initialised as 0)
                        sigma2_d = zeros(self.p,1);                     % the sum of sigma2[d, c] for a given d (must be initialised as 0)
                        % we consider all x[e] -> y[d]
                        for p_id=1:self.p
                            hi = self.his(p_id);
                            li = self.lis(p_id);
                            ki = self.kis(p_id);
                            % calculate coordinates for x[e]
                            e_l = l - 1 - li + 1;
                            if l-1 < li
                                e_l = e_l + self.M;
                            end
                            e_k = mod(k - 1 - ki, self.N) + 1;
                            if self.y_rg.isPulseIdeal()
                                e_h = hi*exp(2j*pi*li/self.M*ki/self.N);
                            elseif self.y_rg.isPulseRecta()
                                e_h = hi*exp(2j*pi*(l-li-1)/self.M*ki/self.N);
                                if l-1 < li
                                    e_h = e_h*exp(-2j*pi*(e_k-1)/self.N);
                                end
                            end
                            % jump if x[e] in PG
                            if pg_num>0 && e_l>=pg_l_beg && e_l<=pg_l_end && e_k>=pg_k_beg && e_k <= pg_k_end
                                continue;
                            end
                            % for the current x[c], we consider all constel points
                            for i2=1:1:self.constel_len
                                mu_d(p_id) = mu_d(p_id) + p_dc(d,p_id,i2) * self.constel(i2);
                                sigma2_d(p_id) = sigma2_d(p_id) + p_dc(d,p_id,i2) * abs(self.constel(i2))^2;
                            end
                            mu_d(p_id) = mu_d(p_id) * e_h;
                            sigma2_d(p_id) = sigma2_d(p_id) * abs(e_h)^2  - abs(mu_d(p_id))^2;
                        end
                        mu_d_sum = sum(mu_d);
                        sigma2_d_sum = sum(sigma2_d)+self.No;
                        % remove x[c] from the sum
                        for p_id=1:self.p
                            mu_dc(d,p_id) = mu_d_sum - mu_d(p_id);
                            sigma2_dc(d,p_id) = sigma2_d_sum - sigma2_d(p_id);
                        end
                    end
                end
                % Update probabilities for each c in p[d, c]
                p_cj = zeros(self.N*self.M, self.constel_len);                     % Pc(aj): the probability that x[c] for each constel point
                for l=1:1:self.M
                    for k=1:1:self.N
                        c = self.N*(l-1)+k;
                        % jump if x[c] in PG area
                        if pg_num>0 && l>=pg_l_beg && l<=pg_l_end && k>=pg_k_beg && k <= pg_k_end
                            continue;
                        end
        
                        pr_ecj_ln = zeros(self.p, self.constel_len);               % ln(Pr(y[e]|x[c]=aj, H))
                        % we consider all y[e] from x[c]
                        e_ls = zeros(self.p,1);                              % y[e] coordinate l
                        e_ks = zeros(self.p,1);                              % y[e] coordinate k
                        for p_id=1:self.p
                            hi = self.his(p_id);
                            li = self.lis(p_id);
                            ki = self.kis(p_id);
                            % y[e] - calculate coordinates 
                            e_l = l - 1 + li + 1;
                            if l + li > self.M
                                e_l = e_l - self.M;
                            end
                            e_k = mod(k - 1 + ki, self.N) + 1;
                            if self.y_rg.isPulseIdeal()
                                e_h = hi*exp(2j*pi*li/self.M*ki/self.N);
                            elseif self.y_rg.isPulseRecta()
                                if l + li <= self.M
                                    e_h = hi*exp(2j*(pi/self.M)*(l-1)*ki/self.N);
                                else
                                    e_h = hi*exp(2j*(pi/self.M)*(l-1-self.M)*ki/self.N)*exp(-2j*pi*(k-1)/self.N);
                                end
                            end
                            % y[e] - record coordinates
                            e_ls(p_id) = e_l;
                            e_ks(p_id) = e_k;
                            % calculate ln(eta(e,c,k)): the probability that x[c] takes a constel point(k) based on y[e]
                            eta_ec_ln = zeros(self.constel_len,1);                 % ln(eta(e,c,k))
                            for i2=1:1:self.constel_len
                                eta_ec_ln(i2) = -(abs(Y_DD(e_k, e_l)- mu_dc(self.N*(e_l-1)+e_k,p_id) - e_h*self.constel(i2))^2)/sigma2_dc(self.N*(e_l-1)+e_k,p_id);
                            end
                            eta_ec_ln = eta_ec_ln - max(eta_ec_ln);     % subtract the maximal exponenet of ln(eta(e,c,k)) to keep mathematical stability
                            % calculate the probability that p[d, c] for the given d
                            pr_ecj_ln(p_id, :) = eta_ec_ln - log(sum(exp(eta_ec_ln)));
                        end
                        % calculate sum(ln(Pr(y[e]|x[c]=aj, H))) for e in J(c)
                        pr_cj_ln = zeros(1, self.constel_len);                     % sum(ln(Pr(y[e]|x[c]=aj, H)))
                        for i2=1:1:self.constel_len
                            pr_cj_ln(i2) = sum(pr_ecj_ln(:,i2));
                        end
                        % calculate p_cj = exp(pr_cj_ln)
                        p_cj_all = exp(pr_cj_ln - max(pr_cj_ln));       % stay mathematical stability
                        p_cj(c, :) = p_cj_all/sum(p_cj_all);            % normalise the sum of probability
                        % calculate p_dc
                        for p_id=1:1:self.p
                            e_l = e_ls(p_id);
                            e_k = e_ks(p_id);
                            % remove ln(Pr(y[d]|x[c]=aj, H)) from sum(ln(Pr(y[e]|x[c]=aj, H)))
                            p_dc_val = pr_cj_ln - pr_ecj_ln(p_id,:);
                            p_dc_val = exp(p_dc_val - max(p_dc_val));
                            p_dc_val = p_dc_val/sum(p_dc_val);
                            p_dc(self.N*(e_l-1)+e_k,p_id,:) = p_dc_val*self.mp_base_delta_fra + (1-self.mp_base_delta_fra)*reshape(p_dc(self.N*(e_l-1)+e_k,p_id,:),1,self.constel_len);
                        end
        
                    end
                end
                % early stop
                conv_rate =  sum(max(p_cj,[],2)>0.99)/(self.N*self.M);
                if conv_rate==1
                    sum_prob_fin = p_cj;
                    break;
                elseif conv_rate > conv_rate_prev
                    conv_rate_prev = conv_rate;
                    sum_prob_fin = p_cj;
                elseif (conv_rate < conv_rate_prev - 0.2) && conv_rate_prev > 0.95
                    break;
                end
            end
            syms = zeros(self.M*self.N - pg_num, 1);
            sym_id = 1;
            for k=1:1:self.N
                for l=1:1:self.M
                    % jump if in PG area
                    if pg_num>0 && l>=pg_l_beg && l<=pg_l_end && k>=pg_k_beg && k <= pg_k_end
                        continue;
                    end
                    [~,pos] = max(sum_prob_fin(self.N*(l-1)+k,:));
                    syms(sym_id) = self.constel(pos);
                    sym_id = sym_id + 1;
                end
            end
        end
    end
end