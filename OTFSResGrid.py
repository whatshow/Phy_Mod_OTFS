import numpy as np
from whatshow_toolbox import MatlabFuncHelper

class OTFSResGrid(MatlabFuncHelper):
    ###########################################################################
    # Constants
    # pulse
    PULSE_NO = 0;
    PULSE_IDEAL = 10;                                   # using ideal pulse to estimate the channel
    PULSE_RECTA = 20;                                   # using rectangular waveform to estimate the channel
    # pilot types
    PILOT_TYPE_EM = 10;                                 # embedded pilots
    PILOT_TYPE_SP = 20;                                 # superimposed pilots
    # Pilot locations
    PILOT_LOC_FLEX = 0;                                 # flexible location
    PILOT_LOC_CENTER = 10;                              # the pilot is put at the center of frame
    PILOT_LOC_ZP = 20;                                  # the pilot is put at the zero padding area
    
    ###########################################################################
    nSubcarNum = 0;                                     # subcarrier number
    nTimeslotNum = 0;                                   # timeslot number
    content = None;
    # zero padding
    zp_len = 0;
    # pulse
    pulse_type = PULSE_NO;
    # pilot
    pilot_type = PILOT_TYPE_EM;
    pilots = [];
    p_len = 0;
    pl1 = 0;                                            # 1st (lowest) pilot location in delay axis
    pk1 = 0;                                            # 1st (lowest) pilot location in Doppler axis
    pl_len = 0;                                         # pilots number along the delay(Doppler) axis
    pk_len = 0;                                         # pilots number along the Doppler(time) axis
    pilot_loc_type = PILOT_LOC_CENTER;                  # pilot attribute
    # guard
    gl_len_ful = False;
    gl_len_neg = 0;
    gl_len_pos = 0;
    gk_len_ful = False;
    gk_len_neg = 0;
    gk_len_pos = 0;
    # channel estimation area in content
    pg_num = 0;
    pg_delay_beg = None;
    pg_delay_end = None;
    pg_doppl_beg = None;
    pg_doppl_end = None;
    # channel estimation area in content
    ce_num = 0;
    ce_delay_beg = None;
    ce_delay_end = None;
    ce_doppl_beg = None;
    ce_doppl_end = None;
    
    ###########################################################################
    # general methods
    '''
    init the resource grid
    @in1:               1st input, a scalar for subcarrier number or the content directly
    @in2:               only if 1st input is scalar, this input is the timeslot number
    @zp_len:            zero padding length
    '''
    def __init__(self, in1, *args, zp_len=0):
        if isinstance(in1, int):
            if len(args) == 0:
                pass
        pass
    
    
    ###########################################################################
    # 