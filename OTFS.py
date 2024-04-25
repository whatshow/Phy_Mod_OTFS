import numpy as np
from numpy import floor, sqrt, exp
from numpy.fft import fft, ifft 
from whatshow_toolbox import MatlabFuncHelper

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
    cp_len = 0;
    
    ###########################################################################
    # General Methods
    '''
    modulate (use fast method by default)
    @rg:        an OTFS resource grid
    @isFast:    DD domain -> TD domain (no X_TF) 
    '''
    def modulate(self, rg, *, isFast=True):
        if not isinstance(rg, 'OTFSResGrid'):
            raise Exception("The input must be an OTFS resource grid.");
        # load RG
        self.nSubcarNum, self.nTimeslotNum = rg.getContentSize();
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
                s_mat = ifft(self.X_TF, axis=-2)*np.sqrt(self.nSubcarNum); # Heisenberg transform
                self.s = self.reshape(s_mat, self.nSubcarNum*self.nTimeslotNum, order='F');
