import numpy as np
from whatshow_toolbox import MatlabFuncHelper
from OTFSResGrid import OTFSResGrid
# constants
pi = np.pi;

class OTFSDetector(MatlabFuncHelper):
    # constants
    # Detect
    DETECT_MP_BASE = 10;                            # base OTFS MP detector proposed by P. Raviteja in 2018
    DETECT_TYPES = [DETECT_MP_BASE];
    #
    constel = None;          # constel
    constel_len = 0;
    M = 0;                   # subcarrier number
    N = 0;                   # timeslot number
    detect_type = None;
    x_num = 0;
    y_rg = None;
    y = None;
    y_num = 0;
    # CSI
    HDD = None;
    p = 0;
    his = [];
    lis = [];
    kis = [];
    No = None;
    # detectors
    # detector - mp base
    mp_base_n_ite = 200;
    mp_base_delta_fra = 0.6;
    
    ###########################################################################
    # General Methods
    '''
    constructor
    @constel:           the constellation (a vector)
    '''
    def __init__(self, constel, *, batch_size=None):
        # constel
        constel = np.asarray(constel).squeeze();
        if constel.ndim != 1:
            raise Exception("The constel must be a vector.");
        else:
            self.constel = constel;
            self.constel_len = len(constel);
        if batch_size is not None:
            self.batch_size = batch_size;
            
    '''
    detect
    @y:             the received signal from the resource grid (after demapping) or just a resource grid 
    @csi_info1:     (1) a vector of path gains  (2) a matrix of HDD
    @csi_info2:     (1) the delay indices       (2) the noise power
    @csi_info3:     (1) the Doppler indices
    @csi_info4:     (1) the noise power
    '''
    def detect(self, y, csi_info1, *args):
        # input
        # input - y
        if isinstance(y, OTFSResGrid):
            self.y_rg = y;
            # nSubcarNum & nTimeslotNum
            [self.M, self.N] = y.getContentSize();
            if self.M <= 0:
                raise Exception("Subcarrier number must be positive.");
            elif self.M != np.floor(self.M):
                raise Exception("Subcarrier number must be an integer.");
            if self.N <= 0:
                raise Exception("Timeslot number must be positive.");
            elif self.N != np.floor(self.N):
                raise Exception("Timeslot number must be an integer.");
            self.y_num = self.M*self.N;
        elif self.isvector(y):
            self.y = y;
            self.y_num = y.shape[-1];
        else:
            raise Exception("y must be a vector from the resource grid (after demapping) or an OTFS resource grid.");
        
        csi_info_extra_len = len(args);
        # input - csi
        csi_info1 = np.asarray(csi_info1);
        if self.isvector(csi_info1):
            # csi - type 1
            self.p = csi_info1.shape[-1];
            self.his = csi_info1;
            if csi_info_extra_len < 2:
                raise Exception("Still need at least two more inputs (the delay indices and the Doppler indices).");
            else:
                self.lis = np.asarray(args[0]);
                self.kis = np.asarray(args[1]);
                if not self.isvector(self.lis) or not self.isvector(self.kis) or self.p != self.lis.shape[-1] or self.p != self.kis.shape[-1]:
                    raise Exception("The CSI must have the same length.");
            if csi_info_extra_len >= 3:
                self.No = np.asarray(args[2]);
        elif self.ismatrix(csi_info1):
            # csi - type 2
            self.HDD = csi_info1;
            y_num_HDD = self.HDD[-2];
            self.x_num = self.HDD[-1];
            if y_num_HDD != self.y_num:
                raise Exception("The HDD row number does not equal to the the received signal number.");
            if csi_info_extra_len >= 1:
                self.No = np.asarray(args[0]);
        else:
            raise Exception("The 1st CSI parameter must be a vector of path gains or a matrix of HDD.");
        # input - csi - No
        if self.No is not None:
            if self.No.ndim != 0:
                raise Exception("The noise power must be a scalar.");
        # detect
        # detect - mp base
        if self.detect_type == self.DETECT_MP_BASE:
            syms = self.detectMPBase();
        return syms;
            
    ###########################################################################
    # Getters & Setters
    '''
    set detector types - MP Base
    @n_ite:           the iteration number (200 by default)
    @delta_fra:       the percentage for taking the values in the current iteration
    '''
    def useMPBase(self, *, n_ite=None, delta_fra=None):
        self.detect_type = self.DETECT_MP_BASE; 
        # MP Base settings
        if n_ite is not None:
            self.mp_base_n_ite = n_ite;
        if delta_fra is not None:
            self.mp_base_delta_fra = delta_fra;
            
    ###########################################################################
    # detectors
    '''
    MP base (proposed by P. Raviteja in 2017) from Emanuele Viterbo Research Group
    '''
    def detectMPBase(self):
        # input check
        if self.No is None or self.No.ndim != 0:
            raise Exception("The noise power must be a scalar.");
        # set data
        conv_rate_prev = -0.1;
        Y_DD = self.y_rg.getContent();
        pg_num, pg_l_beg, pg_l_end, pg_k_beg, pg_k_end = self.y_rg.getAreaPG();   # PG  
        ce_num, ce_l_beg, ce_l_end, ce_k_beg, ce_k_end = self.y_rg.getAreaCE();   # CE
    
        # init all parameters
        # mean & variance (d, c), d=M*(l-1)+k, c=M*(l-li-1)+(k-ki)
        mu_dc = self.zeros(self.N*self.M, self.p).astype(complex);
        sigma2_dc = self.zeros(self.N*self.M, self.p).astype(complex);
        # probability (d, c), d=M*(l-1)+k, c=M*(l-1+li)+(k+ki)
        p_dc = self.ones(self.N*self.M, self.p, self.constel_len)*(1/self.constel_len);
        for ite in range(self.mp_base_n_ite):
            # we update the mean and variance for each d in mu[d, c]
            for l in range(self.M):
                for k in range(self.N):
                    d = self.N*l+k;
                    # jump if y[d] in CE
                    if ce_num>0 and l>=ce_l_beg and l<=ce_l_end and k>=ce_k_beg and k<=ce_k_end:
                        continue;
                    mu_d = self.zeros(self.p).astype(complex);      # the sum of mu[d, c] for a given d (must be initialised as 0)
                    sigma2_d = self.zeros(self.p).astype(complex);  # the sum of sigma2[d, c] for a given d (must be initialised as 0)
                    # we consider all x[e] -> y[d]
                    for p_id in range(self.p):
                        hi = self.his[..., p_id];
                        li = self.lis[..., p_id];
                        ki = self.kis[..., p_id];
                        # calculate coordinates for x[e]
                        e_l = l - li;
                        if l < li:
                            e_l = e_l + self.M;
                        e_k = np.mod(k - ki, self.N);
                        if self.y_rg.isPulseIdeal():
                            e_h = hi*np.exp(2j*pi*li/self.M*ki/self.N);
                        elif self.y_rg.isPulseRecta():
                            e_h = hi*np.exp(2j*pi*(l-li)/self.M*ki/self.N);
                            if l < li:
                                e_h = e_h*np.exp(-2j*pi*e_k/self.N);
                        # jump if x[e] in PG
                        if pg_num>0 and e_l>=pg_l_beg and e_l<=pg_l_end and e_k>=pg_k_beg and e_k <= pg_k_end:
                            continue;
                        # for the current x[c], we consider all constel points
                        for i2 in range(self.constel_len):
                            mu_d[..., p_id] = mu_d[..., p_id] + p_dc[d,p_id,i2]*self.constel[i2];
                            sigma2_d[..., p_id] = sigma2_d[..., p_id] + p_dc[d,p_id,i2]*abs(self.constel[i2])**2;
                        mu_d[..., p_id] = mu_d[..., p_id]*e_h;
                        sigma2_d[..., p_id] = sigma2_d[..., p_id]*abs(e_h)**2  - abs(mu_d[..., p_id])**2;
                    mu_d_sum = self.sum(mu_d);
                    sigma2_d_sum = self.sum(sigma2_d)+self.No;
                    # remove x[c] from the sum
                    for p_id in range(self.p):
                        mu_dc[...,d,p_id] = mu_d_sum - mu_d[..., p_id];
                        sigma2_dc[..., d,p_id] = sigma2_d_sum - sigma2_d[..., p_id];
            # Update probabilities for each c in p[d, c]
            p_cj = self.zeros(self.N*self.M, self.constel_len);     # Pc(aj): the probability that x[c] for each constel point
            # eta_ec_lns = self.zeros(self.M*self.N*self.p, self.constel_len); # debug
            # eta_ec_lns_id = 0;# debug
            for l in range(self.M):
                for k in range(self.N):
                    c = self.N*l+k;
                    # jump if x[c] in PG area
                    if pg_num>0 and l>=pg_l_beg and l<=pg_l_end and k>=pg_k_beg and k <= pg_k_end:
                        continue;
                    pr_ecj_ln = self.zeros(self.p, self.constel_len);    # ln(Pr(y[e]|x[c]=aj, H))
                    # we consider all y[e] from x[c]
                    e_ls = self.zeros(self.p).astype(int);               # y[e] coordinate l
                    e_ks = self.zeros(self.p).astype(int);               # y[e] coordinate k
                    for p_id in range(self.p):
                        hi = self.his[..., p_id];
                        li = self.lis[..., p_id];
                        ki = self.kis[..., p_id];
                        # y[e] - calculate coordinates 
                        e_l = l + li;
                        if l + li >= self.M:
                            e_l = e_l - self.M;
                        e_k = np.mod(k + ki, self.N);
                        if self.y_rg.isPulseIdeal():
                            e_h = hi*np.exp(2j*pi*li/self.M*ki/self.N);
                        elif self.y_rg.isPulseRecta():
                            if l + li < self.M:
                                e_h = hi*np.exp(2j*(pi/self.M)*l*ki/self.N);
                            else:
                                e_h = hi*np.exp(2j*(pi/self.M)*(l-self.M)*ki/self.N)*np.exp(-2j*pi*k/self.N);
                        # y[e] - record coordinates
                        e_ls[..., p_id] = e_l;
                        e_ks[..., p_id] = e_k;
                        # calculate ln(eta(e,c,k)): the probability that x[c] takes a constel point(k) based on y[e]
                        eta_ec_ln = self.zeros(self.constel_len);               # ln(eta(e,c,k))
                        for i2 in range(self.constel_len):
                            eta_ec_ln[...,i2] = np.real(-(abs(Y_DD[..., e_k, e_l]- mu_dc[...,self.N*e_l+e_k,p_id] - e_h*self.constel[i2])**2)/sigma2_dc[self.N*e_l+e_k,p_id]);
                        eta_ec_ln = eta_ec_ln - self.max(eta_ec_ln);            # subtract the maximal exponenet of ln(eta(e,c,k)) to keep mathematical stability
                        # eta_ec_lns[eta_ec_lns_id, :] = eta_ec_ln; # debug
                        # eta_ec_lns_id = eta_ec_lns_id+1;# debug
                        # calculate the probability that p[d, c] for the given d
                        pr_ecj_ln[..., p_id, :] = eta_ec_ln - np.log(self.sum(np.exp(eta_ec_ln)));
                    # calculate sum(ln(Pr(y[e]|x[c]=aj, H))) for e in J(c)
                    pr_cj_ln = self.zeros(1, self.constel_len);                 # sum(ln(Pr(y[e]|x[c]=aj, H)))
                    for i2 in range(self.constel_len):
                        pr_cj_ln[..., i2] = self.sum(pr_ecj_ln[..., :,i2]);
                    # calculate p_cj = exp(pr_cj_ln)
                    p_cj_all = np.exp(pr_cj_ln - self.max(pr_cj_ln));           # stay mathematical stability
                    p_cj[..., c, :] = p_cj_all/self.sum(p_cj_all);              # normalise the sum of probability
                    # calculate p_dc
                    for p_id in range(self.p):
                        e_l = e_ls[p_id];
                        e_k = e_ks[p_id];
                        # remove ln(Pr(y[d]|x[c]=aj, H)) from sum(ln(Pr(y[e]|x[c]=aj, H)))
                        p_dc_val = pr_cj_ln - pr_ecj_ln[..., p_id,:];
                        p_dc_val = np.exp(p_dc_val - self.max(p_dc_val));
                        p_dc_val = p_dc_val/self.sum(p_dc_val);
                        p_dc[..., self.N*e_l+e_k,p_id,:] = p_dc_val*self.mp_base_delta_fra + (1-self.mp_base_delta_fra)*self.reshape(p_dc[..., self.N*e_l+e_k,p_id,:],1,self.constel_len);
            # early stop
            conv_rate =  self.sum(self.max(p_cj)>0.99)/(self.N*self.M);
            if conv_rate == 1:
                sum_prob_fin = p_cj;
                break;
            elif conv_rate > conv_rate_prev:
                conv_rate_prev = conv_rate;
                sum_prob_fin = p_cj;
            elif (conv_rate < conv_rate_prev - 0.2) and conv_rate_prev > 0.95:
                break;
        syms = self.zeros(self.M*self.N - pg_num).astype(complex);
        sym_id = 0;
        for k in range(self.N):
            for l in range(self.M):
                # jump if in PG area
                if pg_num>0 and l>=pg_l_beg and l<=pg_l_end and k>=pg_k_beg and k <= pg_k_end:
                    continue;
                pos = sum_prob_fin[self.N*l+k,:].argmax(axis=-1);
                syms[sym_id] = np.take(self.constel, pos);
                sym_id = sym_id + 1;
        return syms;