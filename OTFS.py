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
    FREQ_SPACING = 15;                          # default frequency spacing is 15kHz
    FC = 3;                                     # default single carrier frequency is 3GHz
    BATCH_SIZE = -1;
    
    # batch
    batch_size = None;                          # the mini batch
    
    # variables
    nSubcarNum = None;                          # subcarrier number
    nTimeslotNum = None;                        # timeslot number
    freq_spacing = None;                        # frequency spacing (kHz), the default is 15kHz
    fc = None;                                  # single carrier frequency (GHz), the default is 3GHz
    X_DD = None;                                # Tx value in the delay Doppler(DD) domain
    X_TF = None;                                # Tx value in the time-frequency(TF) domain
    s = None;                                   # Tx value in the time domain (array)
    H = None;                                   # channel in the time domain
    H_DD = None;                                # Effective channel in the DD domain
    r = None;                                   # Rx value in the time domain (array)
    Y_TF = None;                                # Rx value in the TF domain
    Y_DD = None;                                # Rx value in the DD domain
    taps_num = 0;                               # paths number              
    delay_taps = None;                          # delay index, a row vector
    doppler_taps = None;                        # doppler index (integers or fractional numbers), a row vector
    chan_coef = None;                           # path gain, a row vector
    
    '''
    constructor
    @nSubcarNum:      subcarrier number
    @nTimeslotNum:    timeslot number
    @freq_spacing:    frequency spacing (kHz), the default is 15kHz
    @fc:              single carrier frequency (GHz), the default is 3GHz
    '''
    def __init__(self, nSubcarNum, nTimeslotNum, *, freq_spacing=FREQ_SPACING, fc=FC, batch_size = BATCH_SIZE):
        if not isinstance(nSubcarNum, int):
            raise Exception("The number of subcarriers must be an integer scalar.");
        else:
            self.nSubcarNum = nSubcarNum;
        if not isinstance(nTimeslotNum, int):
            raise Exception("The number of timeslots must be an integer scalar.");
        else:
            self.nTimeslotNum = nTimeslotNum;
        self.freq_spacing = freq_spacing;
        self.fc = fc;
        self.batch_size = batch_size;
    
    '''
    modulate
    @symbols: a vector of symbols to send or a matrix of [Doppler, delay] or [nTimeslotNum ,nSubcarNum]
    '''
    def modulate(self, symbols):
        # input check
        symbols = np.asarray(symbols);
        
    
    def demodulate(self):
        pass
    
    
    def setChannel(self):
        pass
    
    def addChannelPath(self):
        pass
    
    def passChannel(self, noPow):
        pass
    
    # support function
    '''
    Get the channel matrix in Delay Doppler Domain (using the rectangular waveform)
    '''    
    def getChannel(self):
        pass
    '''
    get the signal in the TF domain
    '''
    def getXTF(self):
        pass
    
    '''
    get the signal in the time domain
    '''
    def getS(self):
        pass