import numpy as np
from numpy import floor, sqrt, kron
from numpy.fft import fft, ifft
from whatshow_toolbox import MatlabFuncHelper
from OTFSResGrid import OTFSResGrid
eps = np.finfo(float).eps;

class OTFS(MatlabFuncHelper):
    ###########################################################################
    # Constants
    # Pulse Types
    PULSE_IDEAL = 10;   # ideal pulses (if we use ideal pulses, `ISI_CANCEL_CP` is forced to chosen)
    PULSE_RECTA = 20;   # rectangular pulses
    # cyclic prefix
    CP_ZERO = 10;                                # use zero guards to avoid ISI
    CP_ONE_FRAM = 20;                            # one cp for entire OTFS frame
    CP_ONE_FRAM_SUB = 21;                        # one cp for each OTFS subframe
    
    ###########################################################################
    # system
    fc = 4;                                      # Carrier frequency (GHz)
    fq_sp = 15;                                  # subcarrier spacing (kHz)
    res_k = 0;                                   # resolution - Doppler (kHz)
    res_l = 0;                                   # resolution - delay (ms)
    # RG infomation
    nSubcarNum = None;                           # subcarrier number
    nTimeslotNum = None;                         # timeslot number
    sig_len = 0;                                 # the total signal length
    rg = None
    # OTFS process signal
    X_TF = None                                  # Tx value in the time-frequency(TF) domain
    s  = None;                                   # Tx value in the time domain (array)
    H = None;                                    # channel in the time domain
    r = None;                                    # Rx value in the time domain (array)
    Y_TF = None;                                 # Rx value in the TF domain
    Y_DD = None;
    # channel
    taps_num = 0;                                # paths number           
    chan_coef = None;                            # path gain, a row vector
    delay_taps = None;                           # delay index, a row vector
    doppler_taps = None;                         # doppler index (integers or fractional numbers), a row vector
    # pulse
    pulse_type = PULSE_RECTA;
    # cp
    cp_type = CP_ONE_FRAM;
    cp_len = -0;
    
    ###########################################################################
    # General Methods
    ''''
    constructor
    @batch_size: batch size
    '''
    def __init__(self, *, fc=None, fq_sp=None, batch_size=None):
        if fc is not None:
            self.fc = fc;
        if fq_sp is not None:
            self.fq_sp = fq_sp;
        if batch_size is not None:
            self.batch_size = batch_size;
    
    '''
    modulate (use fast method by default)
    @rg:        an OTFS resource grid
    @isFast:    DD domain -> TD domain (no X_TF) 
    '''
    def modulate(self, rg, *, isFast=True):
        if not isinstance(rg, OTFSResGrid):
            raise Exception("The input must be an OTFS resource grid.");
        # load RG
        self.nSubcarNum, self.nTimeslotNum = rg.getContentSize();
        self.calcRes();
        self.sig_len = self.nSubcarNum*self.nTimeslotNum;
        if rg.isPulseIdeal():
            self.pulse_type = self.PULSE_IDEAL;
        elif rg.isPulseRecta():
            self.pulse_type = self.PULSE_RECTA;
        else:
            raise Exception("Pulse shaping is not given in the resource grid.");
        self.rg = rg.clone();
        X_DD = self.rg.getContent();
        # modulate
        if self.pulse_type == self.PULSE_RECTA:
            if isFast:
                s_mat = ifft(X_DD, axis=-2)*sqrt(self.nTimeslotNum);
                self.s = self.reshape(s_mat, self.nSubcarNum*self.nTimeslotNum);
            else:
                X_FT = fft(ifft(self.X_DD, axis=-2), axis=-1)/sqrt(self.nSubcarNum/self.nTimeslotNum); # ISFFT 
                self.X_TF = np.moveaxis(X_FT, -1, -2); # X_TF is [nSubcarNum, nTimeslotNum]
                s_mat = ifft(self.X_TF, axis=-2)*sqrt(self.nSubcarNum); # Heisenberg transform
                self.s = self.reshape(s_mat, self.nSubcarNum*self.nTimeslotNum, order='F');

    '''
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
    '''
    def setChannel(self, in1, in2, in3, *, force_frac=False, isAWGN=False, isRician=False, isUniform=False):
        in1 = np.asarray(in1);
        in2 = np.asarray(in2);
        in3 = np.asarray(in3);
        # set the channel
        if in1.ndim > 0 and in2.ndim > 0 and in3.ndim > 0:
            # no scalar, giving fixed path
            in1_len = in1.shape[-1];
            in2_len = in2.shape[-1];
            in3_len = in3.shape[-1];
            if not self.isvector(in1) or not self.isvector(in2) or not self.isvector(in3):
                raise Exception("The path gains, delays and dopplers must be 1D vectors.");
            elif in1_len != in2_len or in1_len != in3_len:
                raise Exception("The path gains, delays and dopplers do not have the same length.");
            else:
                self.taps_num = in1_len;
                self.chan_coef = in1;
                self.delay_taps = in2;
                self.doppler_taps = in3;
                self.cp_len = np.max(self.delay_taps).astype(int);
        elif in1.ndim == 0 and in2.ndim == 0 and in3.ndim == 0:
            # scalar, random paths
            p = in1;
            lmax = in2;
            kmax = in3;
            kmax_frac = kmax - floor(kmax);
            kmax = floor(kmax).astype(int);
            if kmax_frac > eps:
                force_frac = True;
            else:
                if force_frac:
                    kmax_frac = 0.5;
            # input check
            if p > (lmax + 1)*(2*kmax+1):
                raise Exception("The path number must be less than lmax*(2*kmax+1) = %d"%(lmax + 1)*(2*kmax+1));
            lmin= 1;
            kmin = -kmax;
            taps_max = (kmax - kmin + 1)*(lmax - lmin + 1);
            # create delay options [lmin, lmin, lmin, lmin+1, lmin+1, lmin+1 ...]
            l_combs = kron(np.arange(1, lmax + 1), np.ones((1, kmax - kmin + 1)));
            # create Doppler options [kmin, kmin+1, kmin+2 ... kmax, kmin ...]
            k_combs = np.tile(np.arange(-kmax, kmax + 1), (1, lmax- lmin + 1));
            # select P paths from all possible paths
            taps_selected_idx = self.shufSelectTopNIdx(taps_max, p);
            # CSI - delay
            self.delay_taps = np.take(l_combs, taps_selected_idx);
            self.delay_taps[..., np.argmin(self.delay_taps, -1)] = 0;
            # CSI - doppler
            self.doppler_taps = k_combs(taps_selected_idx);
            # add fractional Doppler
            if force_frac:
                doppler_taps_k_max_pos_idx = self.doppler_taps == kmax;
                doppler_taps_k_max_neg_idx = self.doppler_taps == -kmax;
                doppler_taps_k_other_idx = abs(self.doppler_taps) != kmax;
                frac_range_max_pos = self.rand(p)*(kmax_frac - kmax + 0.5) - 0.5;
                frac_range_max_neg = self.rand(p)*(kmax - kmax_frac - 0.5) + 0.5;
                frac_range_others = self.rand(p) - 0.5;
                frac_range_all = frac_range_max_pos*doppler_taps_k_max_pos_idx + frac_range_max_neg*doppler_taps_k_max_neg_idx + frac_range_others*doppler_taps_k_other_idx;
                self.doppler_taps = self.doppler_taps + frac_range_all;
            # CSI - others
            self.taps_num = p;
            if isAWGN:
                self.chan_coef = sqrt(1/p)*sqrt(1/2)*(self.randn(p)+1j*self.randn(p));
            elif isRician:
                rician_1st_path = np.asarray([sqrt(1/2/p)]) if self.batch_size == self.BATCH_SIZE_NO else np.asarray([[sqrt(1/2/p)]]*self.batch_size);
                self.chan_coef = np.append(rician_1st_path, sqrt(1/2/p)*(self.randn(p-1)+1j*self.randn(p-1)));
            else:
                self.chan_coef = sqrt(1/2/p)*(self.randn(p)+1j*self.randn(p)); # use Rayleigh fading by default
            self.cp_len = np.max(self.delay_taps).astype(int);
        else:
            raise Exception("The given CSI is not recognised.");
    
    '''
    add a path to the channel (this does not influence other existing paths)
    @hi:      the path gain (linear gain)
    @li:      the delay
    @ki:      the Doppler shift
    '''
    def setChannelExtra(self, hi, li, ki):
        if self.chan_coef is None:
            self.chan_coef = [hi] if self.batch_size ==  self.BATCH_SIZE_NO else hi;
        else:
            self.chan_coef = np.append(self.chan_coef, hi, axis=-1);
        if self.delay_taps is None:
            self.delay_taps = [li] if self.batch_size ==  self.BATCH_SIZE_NO else li;
        else:
            self.delay_taps = np.append(self.delay_taps, li, axis=-1);
        if self.doppler_taps is None:
            self.doppler_taps = [ki] if self.batch_size ==  self.BATCH_SIZE_NO else ki;
        else:
            self.doppler_taps = np.append(self.doppler_taps, ki, axis=-1);
        self.taps_num = self.taps_num + 1;
        self.cp_len = np.max(self.delay_taps);
        
        
    '''
    pass the channel
    @No: noise power (a scalar) or a given noise vector
    '''
    def passChannel(self, No):
        # input check        
        No = self.squeeze(np.asarray(No));
        if No.ndim == 0:
            if No < 0:
                raise Exception("The noise power must be positive.");
        elif self.isvector(No):
            if No.shape[-1] != self.sig_len:
                raise Exception("TThe noise vector length must be %d."%self.sig_len);
        else:
            raise Exception("The noise input must be a scalar for power or a vector for fixed noise.");
        # pass the channel
        s_chan = 0;
        if self.pulse_type == self.PULSE_IDEAL:
            H_DD = self.buildIdealChannel(self.taps_num, self.chan_coef, self.delay_taps, self.doppler_taps);
            s_chan = H_DD @ np.expand_dims(self.rg.getContent(isVector=True), -1);
        elif self.pulse_type == self.PULSE_RECTA:
            s_cp, s_cp_t = self.addCP();
            for tap_id in range(self.taps_num):
                if self.batch_size == self.BATCH_SIZE_NO:
                    hi = self.chan_coef[tap_id]; 
                    li = self.delay_taps[tap_id];
                    ki = self.doppler_taps[tap_id]
                else:
                    hi = np.expand_dims(self.chan_coef[..., tap_id], axis=-1);
                    li =  self.delay_taps[..., tap_id];
                    ki = np.expand_dims(self.doppler_taps[..., tap_id], axis=-1);
                cur_s_tmp = hi*self.circshift(
                    np.append(
                        s_cp*np.exp(2j*np.pi*ki*s_cp_t/self.nSubcarNum/self.nTimeslotNum), 
                        self.zeros(self.cp_len),
                        axis=-1), 
                    li);
                s_chan = s_chan + cur_s_tmp;
            s_chan = self.removeCP(s_chan);
        
        # add noise
        if No.ndim == 0:
            if No > 0:
                noise = sqrt(No/2)*(self.randn(self.sig_len) + 1j*self.randn(self.sig_len));
                s_chan = s_chan + noise;
        elif self.isvector(No):
            s_chan = s_chan + No;
        
        # assign to Rx
        if self.pulse_type == self.PULSE_IDEAL:
            self.Y_DD = self.reshape(s_chan, self.nTimeslotNum, self.nSubcarNum);
        elif self.pulse_type == self.PULSE_RECTA:
            self.r = s_chan;    # to time domain
            return self.r;
        
    '''
    demodulate (use fast method by default)
    @isFast:    TD domain -> DD domain (no Y_TF) 
    '''
    def demodulate(self, *, isFast=True):
        if self.pulse_type == self.PULSE_RECTA:
            r_mat = self.reshape(self.r, self.nSubcarNum, self.nTimeslotNum, order='F');
            if isFast:
                self.Y_DD = fft(np.moveaxis(r_mat, -1, -2), axis=-2)/sqrt(self.nTimeslotNum);
            else:
                # Wigner transform (Y_TF in [nSubcarNum, nTimeslotNum])
                self.Y_TF = fft(r_mat, axis=-2)/sqrt(self.nSubcarNum);
                Y_FT = np.moveaxis(self.Y_TF, -1, -2);
                # SFFT (Y_DD in [Doppler, delay] or [nTimeslotNum ,nSubcarNum])
                self.Y_DD = ifft(fft(Y_FT, axis=-2), axis=-1)/sqrt(self.nTimeslotNum/self.nSubcarNum); 
        self.rg.setContent(self.Y_DD);
        return self.rg;
    
    ###########################################################################
    # Getters & Setters
    '''
    Get the channel matrix in Delay Doppler Domain
    @his:   the channel gains
    @lis:   the channel delays
    @kis:   the channel Dopplers
    @data_only: whether the channel is only for data (by default true). If you want to get the entire H_DD when using pilos and/or guards, you should manullay set it to false.
    '''
    def getChannel(self, *, his=None, lis=None, kis=None, data_only=True):
        # input check & init
        if his is not None:
            his = np.asarray(his);
            lis = np.asarray(lis);
            kis = np.asarray(kis);
            if not self.isvector(his) and not self.isvector(lis) and not self.isvector(kis):
                raise Exception("The input CSI must be vectors.");
            p = his.shape[-1];
            if p != lis.shape[-1] and p != kis.shape[-1]:
                raise Exception("The input CSI (gains, delays and dopplers) must have the same length.");
        else:
            p = self.taps_num;
            his = self.chan_coef;
            lis = self.delay_taps;
            kis = self.doppler_taps;
        # build the channel
        if self.pulse_type == self.PULSE_IDEAL:
            H_DD = self.buildIdealChannel(p, his, lis, kis);
        elif self.pulse_type == self.PULSE_RECTA:
            H_DD = self.buildRectaChannel(p, his, lis, kis);
        # remove the channel for PG & CE
        if data_only:
            H_DD = self.removeNoDAChannel(H_DD);
            
    '''
    get the channel state information
    @sort_by_gain: sort axis
    @sort_by_delay_doppler: sort axes
    @sort_by_doppler_delay: sort axes
    @descend: sort direction
    '''
    def getCSI(self, *, sort_by_gain=False, sort_by_delay_doppler=False, sort_by_doppler_delay=False, descend=False):
        # input check
        if sort_by_gain + sort_by_delay_doppler + sort_by_doppler_delay > 1:
            raise Exception("Cannot sort following over two orders.");
        
        # retrieve CSI
        his = self.chan_coef;
        lis = self.delay_taps;
        kis = self.doppler_taps;
        # TODO: sort
        return his, lis, kis;
    
    '''
    get the signal in the TF domain
    '''
    def getXTF(self):
        return self.X_TF;
    
    '''
    get the signal in the time domain
    '''
    def getXT(self, *, fft_size=None):
        s = None;
        # if fft resolution is lower than subcarrier number, we choose the subcarrier number as the resolution
        if fft_size < self.nSubcarNum:
            s = self.s;
        else:
            # Heisenberg transform
            s_mat = ifft(self.X_TF, n=fft_size, axis=-2)*np.sqrt(self.nSubcarNum);
            # vectorize
            s = self.reshape(s_mat, self.nSubcarNum*self.nTimeslotNum, order='F');
        return s;
    
    ###########################################################################
    # private methods
    '''
    calculate the resolution
    '''
    def calcRes(self):
        self.res_k = self.fq_sp/self.nTimeslotNum;
        self.res_l = 1/self.fq_sp/self.nTimeslotNum;
    
    '''
    shuffle and select top n elements' indices 
    '''
    def shufSelectTopNIdx(self, taps_max, p):
        if self.batch_size == self.BATCH_SIZE_NO:
            taps_idx_chaotic = np.random.permutation(taps_max);
            taps_selected_idx = np.take(taps_idx_chaotic, np.arange(p));
        else:
            taps_selected_idx = np.zeros((self.batch_size, p)).astype(int);
            for batch_id in range(self.batch_size):
                taps_idx_chaotic = np.random.permutation(taps_max);
                taps_selected_idx[batch_id, :] = np.take(taps_idx_chaotic, np.arange(p));
        return taps_selected_idx;
    
    '''
    add cyclic prefix
    '''
    def addCP(self):
        if self.cp_type == self.CP_ZERO:
            s_cp = self.s;
            s_cp_t = self.seq(self.sig_len);
        elif self.cp_type == self.CP_ONE_FRAM:
            s_cp = None;
            if self.batch_size == self.BATCH_SIZE_NO:
                s_cp = np.append(self.s[self.sig_len - self.cp_len:], self.s);
            else:
                s_cp = np.append(self.s[..., self.sig_len - self.cp_len:], self.s, axis=-1);
            s_cp_t = self.seq(-self.cp_len, self.sig_len);
        elif self.cp_type == self.CP_ONE_FRAM_SUB:
            s_mat = self.reshape(self.s, self.nSubcarNum, self.nTimeslotNum, order="F");
            if self.batch_size == self.BATCH_SIZE_NO:
                s_mat = np.append(s_mat[self.nSubcarNum - self.cp_len:, :], s_mat, axis=-2);
            else:
                s_mat = np.append(s_mat[..., self.nSubcarNum - self.cp_len:, :], s_mat, axis=-2);
            s_cp = self.reshape(s_mat, (self.nSubcarNum)*self.nTimeslotNum, order="F");
            s_cp_t = self.repmat1(self.seq(-self.cp_len, self.nSubcarNum), self.nTimeslotNum);
        return s_cp, s_cp_t;
    
    '''
    remove cyclic prefix
    @s_chan: channel output
    '''
    def removeCP(self, s_chan):
        if self.cp_type == self.CP_ZERO:
            s_chan_rm = s_chan;
        elif self.cp_type == self.CP_ONE_FRAM:
            if self.batch_size == self.BATCH_SIZE_NO: 
                s_chan_rm = s_chan[self.cp_len:self.cp_len+self.sig_len];
            else:
                s_chan_rm = s_chan[..., self.cp_len:self.cp_len+self.sig_len];
        elif self.cp_type == self.CP_ONE_FRAM_SUB:
            s_chan_mat = self.reshape(s_chan, -1, self.nTimeslotNum, order="F");
            if self.batch_size == self.BATCH_SIZE_NO: 
                s_chan_rm_mat = s_chan_mat[self.cp_len:self.cp_len+self.nSubcarNum, :];
            else:
                s_chan_rm_mat = s_chan_mat[..., self.cp_len:self.cp_len+self.nSubcarNum, :];
            s_chan_rm = self.reshape(s_chan_rm_mat, self.sig_len, order="F");
        return s_chan_rm;
    
    '''
    remove non-data from DD channel
    @H_DD: the DD channel from
    '''
    def removeNoDAChannel(self, H_DD):
        #TODO: remove zero padding area
        nans = self.nan(self.sig_len);
        pg_num, pg_delay_beg, pg_delay_end, pg_doppl_beg, pg_doppl_end = self.rg.getAreaPG();
        ce_num, ce_delay_beg, ce_delay_end, ce_doppl_beg, ce_doppl_end = self.rg.getAreaCE();
        # mark redundant values - columns (PG area)
        if pg_num > 0:
            for doppl_id in range(pg_doppl_beg, pg_doppl_end+1):
                for delay_id in range(pg_delay_beg, pg_delay_end+1):
                    col_id = doppl_id*self.nSubcarNum + delay_id;
                    H_DD[..., :, col_id] = nans;
        # mark redundant values - rows (CE area)
        if ce_num > 0:
            for doppl_id in range(ce_doppl_beg,ce_doppl_end+1):
                for delay_id in range(ce_delay_beg,ce_delay_end+1):
                    row_id = doppl_id*self.nSubcarNum + delay_id;
                    H_DD[..., row_id, :] = nans;
        # remove
        # remove - columns
        if pg_num > 0:
            col_idx = self.sum(self.isnan(H_DD)) == self.sig_len;
            H_DD = np.delete(H_DD, col_idx, axis=-1);
        # remove - rows
        if ce_num > 0:
            row_idx = self.sum(self.isnan(H_DD), axis=-1) == (self.sig_len - pg_num);
            H_DD = np.delete(H_DD, row_idx, axis=-2);
        return H_DD;
    
    '''
    build the ideal pulse DD channel (callable after modulate)
    @taps_num:  the number of paths
    @his:       the channel gains
    @lis:       the channel delays
    @kis:       the channel dopplers
    '''
    def buildIdealChannel(self, p, his, lis, kis):
        # input check
        if self.pulse_type != self.PULSE_IDEAL:
            raise Exception("Cannot build the ideal pulse DD channel while not using ideal pulse.");
        hw = self.zeros(self.nTimeslotNum, self.nSubcarNum).astype(complex);
        H_DD = self.zeros(self.sig_len, self.sig_len).astype(complex);
        for l in range(self.nSubcarNum):
            for k in range(self.nTimeslotNum):
                    for tap_id in range(p):
                        if self.batch_size == self.BATCH_SIZE_NO:
                            hi = self.chan_coef[tap_id];
                            li = self.delay_taps[tap_id];
                            ki = self.doppler_taps[tap_id];
                        else:
                            hi = np.expand_dims(self.chan_coef[..., tap_id], axis=-1);
                            li = np.expand_dims(self.delay_taps[..., tap_id], axis=-1);
                            ki = np.expand_dims(self.doppler_taps[..., tap_id], axis=-1);
                        hw_add = 1/self.sig_len*hi*np.exp(-2j*np.pi*li*ki/self.sig_len)* \
                                np.expand_dims(np.sum(np.exp(2j*np.pi*(l-li)*self.seq(self.nSubcarNum)/self.nSubcarNum), axis=-1), axis=-1)* \
                                np.expand_dims(np.sum(np.exp(-2j*np.pi*(k-ki)*self.seq(self.nTimeslotNum)/self.nTimeslotNum), axis=-1), axis=-1);
                        if self.batch_size == self.BATCH_SIZE_NO:
                            hw[k, l]= hw[k, l] + hw_add;
                        else:
                            hw[..., k, l]= hw[...,k, l] + self.squeeze(hw_add);
                    if self.batch_size == self.BATCH_SIZE_NO:
                        H_DD = H_DD + hw[k, l]*self.kron(self.circshift(self.eye(self.nTimeslotNum), k), self.circshift(self.eye(self.nSubcarNum), l));
                    else:
                        H_DD = H_DD + np.expand_dims(hw[..., k, l], axis=(-1,-2))*self.kron(self.circshift(self.eye(self.nTimeslotNum), k), self.circshift(self.eye(self.nSubcarNum), l));
        return H_DD;
    
    '''
    build the rectangular pulse DD channel (callable after modulate)
    @taps_num:  the number of paths
    @his:       the channel gains
    @lis:       the channel delays
    @kis:       the channel dopplers
    '''
    def buildRectaChannel(self, p, his, lis, kis):
        # input check
        if self.pulse_type != self.PULSE_RECTA:
            raise Exception("Cannot build the rectangular pulse DD channel while not using rectanular pulse.");
        # build H_DD
        H_DD = self.zeros(self.sig_len, self.sig_len);      # intialize the return channel
        dftmat = self.dftmtx(self.nTimeslotNum);            # DFT matrix
        idftmat = np.conj(dftmat);                          # IDFT matrix
        piMat = self.eye(self.sig_len);                     # permutation matrix (from the delay) -> pi
        # accumulate all paths
        for tap_id in range(p):
            if self.batch_size == self.BATCH_SIZE_NO:
                hi = self.chan_coef[tap_id];
                li = self.delay_taps[tap_id];
                ki = self.doppler_taps[tap_id];
                
            else:
                hi = self.chan_coef[..., tap_id];
                li = self.delay_taps[..., tap_id];
                ki = np.expand_dims(self.doppler_taps[..., tap_id], axis=-1);
            # delay
            piMati = self.circshift(piMat, li);
            # Doppler            
            deltaMat_diag = np.exp(2j*np.pi*ki/(self.sig_len)*self.seq(self.sig_len));
            deltaMati = self.diag(deltaMat_diag);
            # Pi, Qi & Ti
            Pi = self.kron(dftmat, self.eye(self.nSubcarNum)) @ piMati;
            Qi = deltaMati @ self.kron(idftmat, self.eye(self.nSubcarNum));
            Ti = Pi @ Qi;
            # add this path
            if self.batch_size == self.BATCH_SIZE_NO:
                H_DD = H_DD + hi*Ti;
            else:
                H_DD = H_DD + hi.reshape(-1, 1, 1)*Ti;
        return H_DD;