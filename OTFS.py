import numpy as np
from numpy.random import default_rng


## OTFS modulation
#   1. modulation:
#   2. demodulation:
#   3. channel_gen:                         generate the time-variant channel
#   4. channel_out:                         output the time series data after the time-variant channel
#   5. channel_eff_gen:                     generate the effective channel for detection
#   6. Support functions:
#        • get_transmission_symbol_num:     return the number of symbols per transmission
#        • set_transmission_symbols:                 
class OTFS(object):
    # constants
    # Pilot locations
    PILOT_LOC_CENTER = 10;                      # the pilot is put at the center of frame
    PILOT_LOC_DELAY_MOST_CENTER = 20;           # the pilot is put at the most delay area center
    PILOT_LOC_TYPES = [PILOT_LOC_CENTER, PILOT_LOC_DELAY_MOST_CENTER];
    # Detect
    DETECT_MP_BASE = 10;                        # base OTFS MP detector proposed by P. Raviteja in 2018
    DETECT_TYPES = [DETECT_MP_BASE];
    # CSI
    DETECT_CSI_PERFECT = 10;                    # perfect CSI
    DETECT_CSI_CE = 20;                         # CSI from channel estimation (its type is based on pilot type)
    DETECT_CSI_IN = 30;                         # CSI from the input
    DETECT_CSI_TYPES = [DETECT_CSI_PERFECT, DETECT_CSI_CE, DETECT_CSI_IN];
    # Batch
    BATCH_SIZE_NO = None;
    
    # batch
    batch_size = None;                          # the mini batch
    
    # variables
    nSubcarNum = None;                          # subcarrier number
    nTimeslotNum = None;                        # timeslot number
    freq_spacing = 15;                          # frequency spacing (kHz), the default is 15kHz
    fc = 3;                                     # single carrier frequency (GHz), the default is 3GHz
    X_DD = None;                                # Tx value in the delay Doppler(DD) domain
    X_TF = None;                                # Tx value in the time-frequency(TF) domain
    s = None;                                   # Tx value in the time domain (array)
    H = None;                                   # channel in the time domain
    r = None;                                   # Rx value in the time domain (array)
    Y_TF = None;                                # Rx value in the TF domain
    Y_DD = None;                                # Rx value in the DD domain
    y_DD = None;                                # Rx value in the DD domain (vectorized)
    taps_num = 0;                               # paths number              
    delay_taps = None;                          # delay index, a row vector
    doppler_taps = None;                        # doppler index (integers or fractional numbers), a row vector
    chan_coef = None;                           # path gain, a row vector
    # invalid area in X_DD (these parameters will be set after `insertPilotsAndGuards`)
    X_DD_invalid_num = 0;
    X_DD_invalid_delay_beg = None;
    X_DD_invalid_delay_end = None;
    X_DD_invalid_doppl_beg = None;
    X_DD_invalid_doppl_end = None;
    # pilot
    pilots = [];                                # pilots value
    pilots_num_delay = 0;                       # pilots number along the delay(Doppler) axis
    pilots_num_doppler = 0;                     # pilots number along the Doppler(time) axis
    pilot_loc_delay_1st = 0;                    # 1st (lowest) pilot location in delay axis
    pilot_loc_doppl_1st = 0;                    # 1st (lowest) pilot location in Doppler axis
    # pilot location
    pilot_loc_type = PILOT_LOC_CENTER;
    # channel estimation
    ce_delay_beg = None;
    ce_delay_end = None;
    ce_doppl_beg = None;
    ce_doppl_end = None;
    ce_delay_taps = [];                         # estimated delay index, a row vector
    ce_doppler_taps = [];                       # estimated doppler index (integers or fractional numbers), a row vector
    ce_chan_coef = [];                          # estimated path gain, a row vector
    
    ###########################################################################
    # General OTFS Methods
    ###########################################################################
    '''
    constructor
    @nSubcarNum:      subcarrier number
    @nTimeslotNum:    timeslot number
    @freq_spacing:    frequency spacing (kHz), the default is 15kHz
    @fc:              single carrier frequency (GHz), the default is 3GHz
    @pilot_loc_type:  pilot type
    '''
    def __init__(self, nSubcarNum, nTimeslotNum, *, freq_spacing=None, fc=None, pilot_loc_type=None, batch_size = BATCH_SIZE_NO):
        if not isinstance(nSubcarNum, int):
            raise Exception("The number of subcarriers must be an integer scalar.");
        else:
            self.nSubcarNum = nSubcarNum;
        if not isinstance(nTimeslotNum, int):
            raise Exception("The number of timeslots must be an integer scalar.");
        else:
            self.nTimeslotNum = nTimeslotNum;
        # freq_spacing & fc
        if freq_spacing is not None:
            self.freq_spacing = freq_spacing;
        if fc is not None:
            self.fc = fc;
        # pillot
        if pilot_loc_type is not None:
            self.pilot_loc_type = pilot_loc_type;
        # batch_size
        self.batch_size = batch_size;
    
    '''
    modulate
    @symbols: a vector of symbols to send or a matrix of [(batch_size), Doppler, delay] or [(batch_size), nTimeslotNum ,nSubcarNum]
    '''
    def modulate(self, symbols):
        # input check
        symbols = np.asarray(symbols);
        if isvector(symbols, batch_size=self.batch_size):
            if symbols.shape[-1] != self.nTimeslotNum*self.nSubcarNum:
                raise Exception("The transmission symbol number must be %d"%(self.nTimeslotNum*self.nSubcarNum));
        elif ismatrix(symbols, batch_size=self.batch_size):
            if symbols.shape[-2] != self.nTimeslotNum or symbols.shape[-1] !=self.nSubcarNum:
                raise Exception("The transmission symbol must be in the shape of (%d, %d)"%(self.nTimeslotNum, self.nSubcarNum));
        else:
            raise Exception("The transmission symbol must be a vector or a matrix");
        # modulate
        # reshape(rowwise) to [Doppler, delay] or [nTimeslotNum ,nSubcarNum]
        if isvector(symbols, batch_size=self.batch_size):
            self.X_DD = reshape(symbols, (self.nTimeslotNum, self.nSubcarNum), batch_size=self.batch_size);
        else:
            self.X_DD = symbols;
        # ISFFT
        X_FT = np.fft.fft(np.fft.ifft(self.X_DD, axis=-2), axis=-1)/np.sqrt(self.nSubcarNum/self.nTimeslotNum);
        self.X_TF = np.moveaxis(X_FT, -1, -2);
        # Heisenberg transform
        s_mat = np.fft.ifft(self.X_TF, axis=-2)*np.sqrt(self.nSubcarNum);
        # Vectorizations
        self.s = reshape(s_mat, self.nSubcarNum*self.nTimeslotNum, batch_size=self.batch_size, order='F');
    
    '''
    demodulate
    '''
    def demodulate(self):
        r_mat = reshape(self.r, (self.nSubcarNum, self.nTimeslotNum), order='F', batch_size=self.batch_size);
        # Wigner transform (Y_TF in [nSubcarNum, nTimeslotNum])
        self.Y_TF = np.fft.fft(r_mat, axis=-2)/np.sqrt(self.nSubcarNum);
        Y_FT = np.moveaxis(self.Y_TF, -1, -2);
        # SFFT (Y_DD in [Doppler, delay] or [nTimeslotNum ,nSubcarNum])
        self.Y_DD = np.fft.ifft(np.fft.fft(Y_FT, axis=-2), axis=-1)/np.sqrt(self.nTimeslotNum/self.nSubcarNum); 
        # return the DD domain vector (i.e., transpose->vectorization or reshape using 'C' order)
        self.y_DD = reshape(self.Y_DD, self.nSubcarNum*self.nTimeslotNum, batch_size=self.batch_size);
        return self.y_DD;
        
    '''
    set channel
    <set random paths>
    @p: the path number
    @lmax: the maxmimal delay index
    @kmax: the maximal Doppler index (can be fractional)
    <set known paths>
    @delays:    the delays
    @dopplers:  the doppler shifts
    @gains:     the path gains
    @rtnch:     whether return the channel
    '''
    def setChannel(self, *, p=0, lmax=0, kmax=0, delays=[], dopplers=[], gains=[]):
        # input check
        if not isinstance(p, int):
            raise Exception("p must be an integer");
        elif not isinstance(lmax, int):
            raise Exception("lmax must be an integer");
        elif not isinstance(kmax, (int, float)):
            raise Exception("kmax must be an integer or a float number");
        delays = np.asarray(delays);
        if not isinstance(delays, np.ndarray) and (delays.ndim == 1 and self.batch_size == OTFS.BATCH_SIZE_NO or delays.ndim == 2 and self.batch_size != OTFS.BATCH_SIZE_NO):
            raise Exception("delays must be a list");
        dopplers = np.asarray(dopplers);
        if not isinstance(dopplers, np.ndarray) and (dopplers.ndim == 1 and self.batch_size == OTFS.BATCH_SIZE_NO or dopplers.ndim == 2 and self.batch_size != OTFS.BATCH_SIZE_NO):
            raise Exception("dopplers must be a list");
        gains = np.asarray(gains);
        if not isinstance(gains, np.ndarray) and (gains.ndim == 1 and self.batch_size == OTFS.BATCH_SIZE_NO or gains.ndim == 2 and self.batch_size != OTFS.BATCH_SIZE_NO):
            raise Exception("gains must be a list");
        if delays.size != dopplers.size or delays.size != gains.size:
            raise Exception("The delays, dopplers and gains do not have the same length.");
        
        # reset paramets if using fractional Doppler
        is_fractional_doppler = False;
        kmax_frac = 0.0;
        if kmax != np.floor(kmax):
            is_fractional_doppler = True;
            kmax_frac = kmax;
            kmax = np.floor(kmax_frac).astype(int);
        
        # set the channel
        if delays.size != 0 and dopplers.size != 0 and gains.size != 0:
            if delays.shape[-1] != dopplers.shape[-1] and delays.shape[-1] != gains.shape[-1]:
                raise Exception("The delays, dopplers and gains do not have the same length.");
            else:
                self.delay_taps = delays;
                self.doppler_taps = dopplers;
                self.chan_coef = gains;
                self.taps_num = delays.shape[-1];
        elif p > 0 and lmax >= 1 and (kmax > 0 or kmax_frac > 0):
            # create delay options [lmin, lmin, lmin, lmin+1, lmin+1, lmin+1 ...]
            l = np.arange(1, lmax + 1);
            l_combs = kron(l, np.ones((1, kmax*2 + 1))).astype(int);             # delays are all integers
            # create Doppler options [kmin, kmin+1, kmin+2 ... kmax, kmin ...]
            k = np.arange(-kmax, kmax + 1);
            k_combs = np.tile(k, (1, lmax + 1));
            # We select P paths from all possible paths; that is, we do the randperm(taps_max) and we choose the first P items
            taps_max = (kmax*2 + 1)*lmax;
            taps_selected_idx = Channl_SelectRandomPathIdx(taps_max, p, batch_size=self.batch_size);
            # set channel info
            self.delay_taps = np.take(l_combs, taps_selected_idx);
            # the 1st minimal delay is 0
            self.delay_taps[..., np.argmin(self.delay_taps, -1)] = 0;
            self.doppler_taps = np.take(k_combs, taps_selected_idx);
            # append fractional Doppler
            if is_fractional_doppler:
                doppler_taps_k_max_pos_idx = self.doppler_taps == kmax;
                doppler_taps_k_max_neg_idx = self.doppler_taps == -kmax;
                doppler_taps_k_other_idx = abs(self.doppler_taps) != kmax;
                frac_range_max_pos = self.rand(p)*(kmax_frac - kmax + 0.5) - 0.5;
                frac_range_max_neg = self.rand(p)*(kmax - kmax_frac - 0.5) + 0.5;
                frac_range_others = self.rand(p) - 0.5;
                frac_range_all = frac_range_max_pos*doppler_taps_k_max_pos_idx + frac_range_max_neg*doppler_taps_k_max_neg_idx + frac_range_others*doppler_taps_k_other_idx;
                self.doppler_taps = self.doppler_taps + frac_range_all;
            self.chan_coef = np.sqrt(1/p)*np.sqrt(1/2)*self.randn(p) + 1j*self.randn(p);
            self.taps_num = p;
        else:
            raise Exception("Channel Infomation is not recognised.");
    
    '''    
    add a path to the channel (this does not influence other existing paths)
    @hi:      the path gain (linear gain)
    @li:      the delay
    @ki:      the Doppler shift
    '''
    def addChannelPath(self, hi, li, ki):
        if self.chan_coef is None:
            self.chan_coef = [hi] if self.batch_size ==  OTFS.BATCH_SIZE_NO else hi;
        else:
            self.chan_coef = np.append(self.chan_coef, hi, axis=-1);
        if self.delay_taps is None:
            self.delay_taps = [li] if self.batch_size ==  OTFS.BATCH_SIZE_NO else li;
        else:
            self.delay_taps = np.append(self.delay_taps, li, axis=-1);
        if self.doppler_taps is None:
            self.doppler_taps = [ki] if self.batch_size ==  OTFS.BATCH_SIZE_NO else ki;
        else:
            self.doppler_taps = np.append(self.doppler_taps, ki, axis=-1);
        self.taps_num = self.taps_num + 1;
    
    '''
    pass the channel
    @No: noise power (a scalar) or a given noise vector
    '''
    def passChannel(self, No):
        No = self.squeeze(np.asarray(No));
        if No.ndim == 0:
            if No < 0:
                raise Exception("The noise power must be positive.");
        elif self.isvector(No):
            if No.shape[-1] != self.nTimeslotNum*self.nSubcarNum:
                raise Exception("TThe noise vector length must be %d."%(self.nSubcarNum*self.nTimeslotNum));
        else:
            raise Exception("The noise input must be a scalar for power or a vector for fixed noise.");
        # add CP
        cp_len = np.max(self.delay_taps);
        s_cp = Channel_AddCP(self.s, cp_len, batch_size=self.batch_size);
        # pass the channel
        s_chan = self.zeros(s_cp.shape[-1] + cp_len);
        for tap_id in range(self.taps_num):
            hi = self.chan_coef[tap_id] if self.batch_size == OTFS.BATCH_SIZE_NO else np.expand_dims(self.chan_coef[..., tap_id], axis=-1);
            li = self.delay_taps[tap_id] if self.batch_size == OTFS.BATCH_SIZE_NO else self.delay_taps[..., tap_id];
            ki = self.doppler_taps[tap_id] if self.batch_size == OTFS.BATCH_SIZE_NO else np.expand_dims(self.doppler_taps[..., tap_id], axis=-1);
            cur_s_tmp = hi*circshift(
                np.append(
                    np.multiply(s_cp, np.exp(1j*2*np.pi/self.nSubcarNum*seq(-cp_len, -cp_len+s_cp.shape[-1]-1, order='F', batch_size=self.batch_size)*ki/self.nTimeslotNum)), 
                    self.zeros(cp_len),
                    axis=-1), 
                li,
                batch_size=self.batch_size);
            s_chan = s_chan + cur_s_tmp;
        self.r = s_chan;
        # remove CP
        self.r = self.r[cp_len:cp_len+(self.nTimeslotNum*self.nSubcarNum)] if self.batch_size == OTFS.BATCH_SIZE_NO else self.r[:, cp_len:cp_len+(self.nTimeslotNum*self.nSubcarNum)];
        # add noise
        if No.ndim == 0:
            if No > 0:
                noise = np.sqrt(No/2)*(self.randn(self.nTimeslotNum*self.nSubcarNum) + 1j*self.randn(self.nTimeslotNum*self.nSubcarNum));
                self.r = self.r + noise;
        elif self.isvector(No):
            self.r = self.r + No;
        # return
        return self.r;
    
    ###########################################################################
    # Channel estimation
    ###########################################################################
    
    ###########################################################################
    # OTFS Detectors
    ###########################################################################
    
    ###########################################################################
    # support function
    ###########################################################################
    '''
    Get the channel matrix in Delay Doppler Domain (using the rectangular waveform)
    '''    
    def getChannel(self):
        # intialize the return channel
        H_DD = self.zeros(self.nTimeslotNum*self.nSubcarNum, self.nTimeslotNum*self.nSubcarNum);
        # DFT & IDFT matrix
        dftmat = dftmtx(self.nTimeslotNum, batch_size=self.batch_size);
        idftmat = np.conj(dftmat);
        # permutation matrix (from the delay) -> pi
        piMat = eye(self.nTimeslotNum*self.nSubcarNum, batch_size=self.batch_size);
        #accumulate all paths
        for tap_id in range(self.taps_num):
            li = self.delay_taps[tap_id] if self.batch_size == OTFS.BATCH_SIZE_NO else self.delay_taps[..., tap_id];
            ki = self.doppler_taps[tap_id] if self.batch_size == OTFS.BATCH_SIZE_NO else self.doppler_taps[..., tap_id];
            hi = self.chan_coef[tap_id] if self.batch_size == OTFS.BATCH_SIZE_NO else self.chan_coef[..., tap_id];
            # delay
            piMati = circshift(piMat, li, batch_size=self.batch_size);
            # Doppler
            deltaMat_diag_idx = seq(self.nTimeslotNum*self.nSubcarNum, batch_size=self.batch_size);
            ki_val = ki;
            if self.batch_size != OTFS.BATCH_SIZE_NO:
                ki_val = np.expand_dims(ki_val, axis=1); # transfer its shape -> (batch_size, 1)
            deltaMat_diag = np.exp(1j*2*np.pi*ki_val/(self.nTimeslotNum*self.nSubcarNum)*deltaMat_diag_idx);
            deltaMati = diag(deltaMat_diag, batch_size=self.batch_size);
            # Pi
            Pi = np.matmul(kron(dftmat, eye(self.nSubcarNum, batch_size=self.batch_size), batch_size=self.batch_size), piMati);
            # Qi
            Qi = np.matmul(deltaMati, kron(idftmat, eye(self.nSubcarNum, batch_size=self.batch_size), batch_size=self.batch_size));
            # generate Ti
            Ti = np.matmul(Pi, Qi);
            # add this path
            if self.batch_size == OTFS.BATCH_SIZE_NO:
                H_DD = H_DD + hi*Ti;
            else:
                H_DD = H_DD + hi.reshape(-1, 1, 1)*Ti;
        # return
        return H_DD;
    
    '''
    get the channel delays
    '''
    def getChannelDelays(self):
        return self.delay_taps;
    
    '''
    get the channel dopplers
    '''
    def getChannelDopplers(self):
        return self.doppler_taps;
    
    '''
    get the channel gains
    '''
    def getChannelGains(self):
        return self.chan_coef;
    
    '''
    get the signal in the Delay Time domain [delay, time]
    '''
    def getX2DT(self):
        return np.moveaxis(np.fft.ifft(self.X_DD), -1, -2);
    
    '''
    get the signal in the TF domain
    '''
    def getX2TF(self):
        return self.X_TF;
    
    '''
    get the signal in the time domain
    '''
    def getX2T(self, *, fft_size=None):
        s = None;
        # if fft resolution is lower than subcarrier number, we choose the subcarrier number as the resolution
        if fft_size < self.nSubcarNum:
            s = self.s;
        else:
            # Heisenberg transform
            s_mat = np.fft.ifft(self.X_TF, n=fft_size, axis=-2)*np.sqrt(self.nSubcarNum);
            # vectorize
            s = reshape(s_mat, self.nSubcarNum*self.nTimeslotNum, batch_size=self.batch_size, order='F');
        return s;
    
    '''
    get the received signal in the delay Doppler domain
    '''
    def getYDD(self):
        return self.Y_DD;
    
    ##########################################################################
    # Functions uniform with non-batch and batch
    ##########################################################################
    '''
    check input is a vector like [(batch_size), n],  [(batch_size), n, 1] or [(batch_size), 1, n] 
    '''
    def isvector(self, mat):
        mat = np.asarray(mat);
        if self.batch_size is self.BATCH_SIZE_NO:
            return mat.ndim == 1 or mat.ndim == 2 and (mat.shape[-2] == 1 or mat.shape[-1] == 1);
        else:
            if mat.shape[0] != self.batch_size:
                raise Exception("The input does not has the required batch size.");
            else:
                return mat.ndim == 2 or mat.ndim == 3 and (mat.shape[-2] == 1 or mat.shape[-1] == 1);
    
    '''
    check input is a matrix like [(batch_size), n. m]
    '''
    def ismatrix(self, mat):
        mat = np.asarray(mat);
        if self.batch_size is self.BATCH_SIZE_NO:
            return mat.ndim == 2 and mat.shape[-2] > 1 and mat.shape[-1] > 1;
        else:
            if mat.shape[0] != self.batch_size:
                raise Exception("The input does not has the required batch size.");
            else:
                return mat.ndim == 3 and mat.shape[-2] > 1 and mat.shape[-1] > 1;
            
    '''
    generate a matrix of all zeros
    @order: 'C': this function only create given dimensions; 'F': create the dimensions as matlab (2D at least)
    '''
    def zeros(self, nrow, *args, order='C'):
        out = None;
        if order == 'F':
            ncol = args[0] if len(args) >= 1 else nrow;
            out = np.zeros((nrow, ncol)) if self.batch_size == OTFS.BATCH_SIZE_NO else np.zeros((self.batch_size, nrow, ncol));
        elif order == 'C':
            zeros_shape = list(args);
            zeros_shape.insert(0, nrow);
            if self.batch_size != OTFS.BATCH_SIZE_NO:
                zeros_shape.insert(0, self.batch_size);
            out = np.zeros(zeros_shape);
        return out;
    
    '''
    generate random values from [0, 1) following a uniform distribution
    @args: d0, d1, ..., dn (multiple inputs)
    '''
    def rand(self, *args):
        if self.batch_size is OTFS.BATCH_SIZE_NO:
            return np.random.rand(*args);
        else:
            return np.random.rand(self.batch_size, *args);
        
    '''
    generate random values following the standard normal distribution 
    @args: d0, d1, ..., dn (multiple inputs)
    '''
    def randn(self, *args):
        if self.batch_size is OTFS.BATCH_SIZE_NO:
            return np.random.randn(*args);
        else:
            return np.random.randn(self.batch_size, *args);
        
    '''
    squeeze redundant dimension except the batch size
    '''
    def squeeze(self, mat):
        out = mat.copy();
        out = np.squeeze(out);
        if self.batch_size == 1 and mat.ndim > 0:
            out = np.expand_dims(out, 0);
        return out;

##############################################################################
# Support Functions (only used in this model, please don't import)
# Dedicated for 2D & 3D
##############################################################################
BATCH_SIZE_NO = OTFS.BATCH_SIZE_NO;
'''
check the input is a vector or not
'''
def isvector(mat, *, batch_size=BATCH_SIZE_NO):
    mat = np.asarray(mat);
    if batch_size is BATCH_SIZE_NO:
        return mat.ndim == 1;
    else:
        return mat.ndim == 2;

'''
check the input is a matrix or not
'''
def ismatrix(mat, *, batch_size=BATCH_SIZE_NO):
    mat = np.asarray(mat);
    if batch_size is BATCH_SIZE_NO:
        return mat.ndim == 2;
    else:
        return mat.ndim == 3;

'''
generate a sequence
'''
def seq(*args, batch_size=BATCH_SIZE_NO, order='C'):
    out = None;
    seq_start = 0 if order == 'C' else 1;    
    seq_step = 1;
    seq_end = None;
    if len(args) == 1:
        seq_end = args[0] if order == 'C' else args[0] + 1;
    elif len(args) == 2:
        seq_start = args[0];
        seq_end = args[1] if order == 'C' else args[1] + 1;
    elif len(args) == 2:
        seq_start = args[0];
        seq_step = args[1];
        seq_end = args[2] if order == 'C' else args[2] + 1;
    out = np.arange(seq_start, seq_end, seq_step);
    if batch_size is not BATCH_SIZE_NO:
        out = np.tile(out,(batch_size, 1));
    return out;

'''
generate identity matrix
'''
def eye(size, *, batch_size=BATCH_SIZE_NO):
    out = np.eye(size);
    if batch_size is not BATCH_SIZE_NO:
        out = np.tile(out,(batch_size, 1, 1));
    return out;

'''
reshape
'''
def reshape(mat, shape, *, order='C', batch_size=BATCH_SIZE_NO):
    if isinstance(shape, int):
        shape = [shape];
    else:
        shape = list(shape);
    if batch_size is not BATCH_SIZE_NO:
        shape.insert(0, batch_size);
    return np.reshape(mat, shape, order=order);


'''
generate the DFT matrix
'''
def dftmtx(size, *, batch_size=BATCH_SIZE_NO):
    ident_mat = eye(size, batch_size=batch_size);
    out = np.fft.fft(ident_mat)/np.sqrt(size);
    return out;

'''
circular shift (1st index except for the batch size)
@mat: [(batch_size), dim1, dim2, ...]
@step: an 1D array when using batches or a scalar when not
'''
def circshift(mat, step, *, batch_size=BATCH_SIZE_NO):
    out = None;
    if batch_size is BATCH_SIZE_NO:
        out = np.roll(mat, step, axis=0);
    else:
        if mat.ndim == 2:
            out = np.zeros((batch_size, mat.shape[-1]), dtype=mat.dtype);
        if mat.ndim == 3:
            out = np.zeros((batch_size, mat.shape[-2], mat.shape[-1]), dtype=mat.dtype);
        for batch_id in range(batch_size):
            out[batch_id, ...] = np.roll(mat[batch_id, ...], step[batch_id], axis=0);
    return out;

'''
generate a matrix based on its diag
'''
def diag(diag_vec, *, batch_size=BATCH_SIZE_NO):
    out = None;
    diag_vec_len = diag_vec.shape[-1];
    if batch_size is BATCH_SIZE_NO:
        out = np.diag(diag_vec);
    else:
        # np.zeros only take real numbers by default, here we need to put complex value into it
        out = np.zeros((batch_size, diag_vec_len, diag_vec_len), dtype=diag_vec.dtype);
        for batch_id in range(batch_size):
            out[batch_id, ...] = np.diag(diag_vec[batch_id, ...]);
    return out;

'''
Kronecker product
@a: a matrix like [(batch_size), i, j]
@b: a matrix like [(batch_size), k, l]
'''
def kron(a, b, *, batch_size=BATCH_SIZE_NO):
    out = None;
    if batch_size is BATCH_SIZE_NO:
        out = np.kron(a, b);
    else:    
        sh=a.shape[-1]*b.shape[-1];
        out = np.einsum('aij,akl->aikjl', a, b).reshape(batch_size, sh, sh);
    return out;

##############################################################################
# Cross Dimension Functions (only used in this model, please don't import)
# Dedicated for 2D & 3D
##############################################################################
'''
generate a sequence from 0 to size-1 in a random sequence
'''
def Channl_SelectRandomPathIdx(taps_max, p, *, batch_size=BATCH_SIZE_NO):
    if batch_size == BATCH_SIZE_NO:
        taps_idx_chaotic = np.random.permutation(taps_max);
        taps_selected_idx = np.take(taps_idx_chaotic, np.arange(p));
    else:
        taps_selected_idx = np.zeros((batch_size, p)).astype(int);
        for batch_id in range(batch_size):
            taps_idx_chaotic = np.random.permutation(taps_max);
            taps_selected_idx[batch_id, :] = np.take(taps_idx_chaotic, np.arange(p));
    return taps_selected_idx;

'''
add cyclic prefix
'''
def Channel_AddCP(s, cp_len, *, batch_size=BATCH_SIZE_NO):
    s_cp = None;
    syms_len = s.shape[-1];
    if batch_size == BATCH_SIZE_NO:
        s_cp = np.append(s[syms_len - cp_len:], s);
    else:
        s_cp = np.append(s[..., syms_len - cp_len:], s, axis=-1);
    return s_cp;