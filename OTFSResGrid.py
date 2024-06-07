import numpy as np
from numpy import pi
from numpy import floor, sqrt, exp
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
    pilots = np.array([]);
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
    @batch_size:        the batch_size
    '''
    def __init__(self, in1, *args, zp_len=0, batch_size=None):
        # optional inputs
        self.zp_len = zp_len;
        if batch_size is not None:
            self.batch_size = batch_size;
        # inputs
        in1 = self.squeeze(in1);
        if in1.ndim == 0:
            if len(args) < 1:
                raise Exception("The timeslot number is not given.");
            else:
                in2 = self.squeeze(args[0]);
                if floor(in1) != in1:
                    raise Exception("The subcarier number can only be an integer.");
                if in2.ndim > 0:
                    raise Exception("The timeslot number is not a scalar.");
                elif floor(in2) != in2:
                    raise Exception("The timeslot number can only be an integer.");
                else:
                    self.nSubcarNum = in1.astype(int);
                    self.nTimeslotNum = in2.astype(int);
                    self.content = self.zeros(self.nTimeslotNum, self.nSubcarNum).astype(np.complex_);
        else:
            if self.batch_size == self.BATCH_SIZE_NO and in1.ndim < 2 or self.batch_size != self.BATCH_SIZE_NO and in1.ndim < 3:
                in1 = np.expand_dims(in1, axis=-2);
            self.content = in1;
            self.nTimeslotNum = in1.shape[-2];
            self.nSubcarNum = in1.shape[-1];
    
    '''
    pulse settings
    '''
    def setPulse2Ideal(self):
        self.pulse_type = self.PULSE_IDEAL;
    def setPulse2Recta(self):
        self.pulse_type = self.PULSE_RECTA;
        
    '''
    pilot type setting
    '''
    def setPilot2Embed(self):
        self.pilot_type = self.PILOT_TYPE_EM;
    def setPilot2SuperImposed(self):
        self.pilot_type = self.PILOT_TYPE_SP;
        
    '''
    pilot position setting
    @pl_len:    pilot length on the delay
    @pk_len:    pilot length on the doppler
    @pl1:       pilot location on the delay
    @pk1:       pilot location on the doppler
    '''
    def setPilot2Center(self, pl_len, pk_len):
        self.pilot_loc_type = self.PILOT_LOC_CENTER;
        self.pl_len = pl_len;
        self.pk_len = pk_len;
        self.pl1 = floor((self.nSubcarNum - self.pl_len)/2).astype(int);
        self.pk1 = floor((self.nTimeslotNum - self.pk_len)/2).astype(int);
        self.validP();
    def setPilot2ZP(self, pl_len, pk_len):
        if self.zp_len == 0:
            raise Exception("Zero Padding is not used in this resource grid.");
        self.pilot_loc_type = self.PILOT_LOC_ZP;
        self.pl_len = pl_len;
        self.pk_len = pk_len;
        self.pl1 = floor((self.nSubcarNum - self.pl_len)/2).astype(int);
        self.pk1 = self.nTimeslotNum - self.zp_len + floor((self.zp_len - self.pk_len)/2).astype(int);
        self.validP();
    def setPilot2Flex(self, pl_len, pk_len, pl1, pk1):
        self.pilot_loc_type = self.PILOT_LOC_FLEX;
        self.pl_len = pl_len;
        self.pk_len = pk_len;
        self.pl1 = pl1;
        self.pk1 = pk1;
        self.validP();
        
    '''
    set guards
    @in1: (1,2) negative guard on the delay (3) negative guard on the Doppler
    @in2: (1,2) positive guard on the delay (3) positive guard on the Doppler
    @in3: (1) negative guard on the Doppler
    @in4: (1) positive guard on the Doppler
    @guard_delay_full:  full guard on delay (if set true, ignore the number setting)
    @guard_doppl_full:  full guard on Doppler (if set true, ignore the number setting)
    '''
    def setGuard(self, *args, guard_delay_full=False, guard_doppl_full=False):
        # pilot check
        if self.pk_len <= 0 or self.pl_len <= 0:
            raise Exception("Pilots must be set before setting guards.");
        # take inputs
        args_len = len(args);
        ins = np.asarray([0,0,0,0]);                                  # guard lengths on 4 directions
        ins_len = args_len;
        for arg_id in range(args_len):
            ins[arg_id] = args[arg_id];
        self.gl_len_ful = guard_delay_full;
        self.gk_len_ful = guard_doppl_full;
        # input check - full guard
        if self.gl_len_ful:
            if ins_len == 2:    # the two inputs are Doppler settings
                ins[2:4] = ins[0:2];
            ins[0] = floor((self.nSubcarNum - self.pl_len)/2);
            ins[1] = self.nSubcarNum - self.pl_len - ins[0];
        if self.gk_len_ful:
            ins[2] = floor((self.nTimeslotNum - self.pk_len)/2);
            ins[3] = self.nTimeslotNum - self.pk_len - ins[2];
        # input check - guard - integers only
        if ins.dtype != np.intc and ins.dtype != np.uintc and ins.dtype != np.int_ :
            raise Exception("Guard number along the delay/Doppler axis must be integers.");
        # input check - guard - no negative
        if ins[0] < 0 or ins[1] < 0 or  ins[2] < 0 or ins[3] < 0:
            raise Exception("Guard number along the delay/Doppler axis must be non-negative.");       
        # take inputs
        self.gl_len_neg = ins[0];
        self.gl_len_pos = ins[1];
        self.gk_len_neg = ins[2];
        self.gk_len_pos = ins[3];
        
    '''
    map
    @symbols:       OTFS symbols
    @pilots:        a vector of your pilots (if given `pilots_pow` won't be used)
    @pilots_pow:    pilot power to generate random pilots
    '''
    def map(self, symbols, *, pilots=[], pilots_pow=None):
        self.pilots = np.asarray(pilots);
        # calculate PG & CE area
        self.calcAreaPGCE();
        # insert
        self.insertDA(symbols);
        self.insertP(pilots_pow);
        
    '''
    set the channel estimate area (CE area is set in `map`, if you want to use your own area, call this)
    @ce_l_beg: CE delay beginning
    @ce_l_end: CE delay ending
    @ce_k_beg: CE Doppler beginning
    @ce_k_end: CE Doppler ending
    '''
    def setAreaCE(self, ce_delay_beg, ce_delay_end, ce_doppl_beg, ce_doppl_end):
        if ce_delay_beg < 0 or ce_delay_beg >= self.nSubcarNum:
            raise Exception("CE delay beginning overflows.");
        elif ce_delay_end < 0 or ce_delay_end >= self.nSubcarNum:
            raise Exception("CE delay ending overflows.");
        elif ce_delay_beg > ce_delay_end:
            raise Exception("CE delay beginning is after CE delay ending.");
        if ce_doppl_beg < 0 or ce_doppl_beg >= self.nTimeslotNum:
            raise Exception("CE Doppler beginning overflows.");
        elif ce_doppl_end < 0 or ce_doppl_end >= self.nTimeslotNum:
            raise Exception("CE Doppler ending overflows.");
        elif ce_doppl_beg > ce_doppl_end:
            raise Exception("CE Doppler beginning is after CE delay ending.");
        self.ce_delay_beg = ce_delay_beg;
        self.ce_delay_end = ce_delay_end;
        self.ce_doppl_beg = ce_doppl_beg;
        self.ce_doppl_end = ce_doppl_end;
        self.ce_num = (ce_delay_end - ce_delay_beg + 1)*(ce_doppl_end - ce_doppl_beg + 1);
    
    '''
    demap
    @isData:        whether give the data
    @isDataVec:     whether the data is vectorized (default: true)
    @isCE:          whether give the channel estimation result
    @threshold:     the threshold to estimate the channel
    '''
    def demap(self, *, isData=True, isDataVec=True, isCE=True, threshold=0):
        # input check          
        # input check - threshold
        if threshold < 0:
            raise Exception("The threshould must be non-negative.");
        # input check - pulse_type
        if self.pulse_type == self.PULSE_NO:
            raise Exception("The pulse type has to be set before demapping.");

        # y
        y = None;
        if isData:
            if self.pilot_type == self.PILOT_TYPE_EM:
                if isDataVec:
                    y = self.getContentNoCE();
                else:
                    y = self.getContentZeroCE();
            elif self.pilot_type == self.PILOT_TYPE_SP:
                    y = self.getContent(isVector=isDataVec);
        his = None;
        lis = None;
        kis = None;
        # need to estimate the channel & have pilots
        if isCE and self.isPG():
            if self.p_len == 1:
                his, lis, kis = self.estimateChannel4SingPilot(threshold);
        return y, his, lis, kis;
    
    ###########################################################################
    # issers, getters, setters
    '''
    clone
    '''
    def clone(self):
        rg = OTFSResGrid(self.content);
        rg.batch_size = self.batch_size;
        rg.zp_len = self.zp_len;
        rg.pulse_type = self.pulse_type;
        rg.pilots = self.pilots;
        rg.p_len = self.p_len;
        rg.pl_len = self.pl_len;
        rg.pk_len = self.pk_len;
        rg.pl1 = self.pl1;
        rg.pk1 = self.pk1;
        rg.pilot_loc_type = self.pilot_loc_type;
        rg.pg_num = self.pg_num;
        rg.pg_delay_beg = self.pg_delay_beg;
        rg.pg_delay_end = self.pg_delay_end;
        rg.pg_doppl_beg = self.pg_doppl_beg;
        rg.pg_doppl_end = self.pg_doppl_end;
        rg.ce_num = self.ce_num;
        rg.ce_delay_beg = self.ce_delay_beg;
        rg.ce_delay_end = self.ce_delay_end;
        rg.ce_doppl_beg = self.ce_doppl_beg;
        rg.ce_doppl_end = self.ce_doppl_end;
        return rg;

    '''
    return zero padding length
    '''
    def isZP(self):
         return self.zp_len;

    '''
    check whether use pilots & guards
    '''
    def isPG(self):
        # check whether pilots is assigned or not
        is_pg = self.pilots.shape[-1] > 0;
        # check whether pilot and guard area is calculated
        is_pg = is_pg and self.pg_num > 0;
        # check whether CE area is calulated
        is_pg = is_pg and self.ce_num > 0;
        return is_pg;

    '''
    check whether the current position is in channel estimation area
    @pos_doppl:     the position on the Doppler axis. Not given means the position is for a Doppler-delay vector
    @pos_delay:     the position on the delay axis for matrix or th position on the Doppler-delay axis for the vector
    '''
    def isInAreaPG(self, pos_doppl, *args):
        return self.isInArea(10, pos_doppl, *args);
    def isInAreaZP(self, pos_doppl, *args):
        return self.isInArea(11, pos_doppl, *args);
    def isInAreaDA(self, pos_doppl, *args):
        return self.isInArea(12, pos_doppl, *args);
    def isInAreaCE(self, pos_doppl, *args):
        return self.isInArea(20, pos_doppl, *args);

    '''
    get the area of pilots and guards
    '''
    def getAreaPG(self):
        return self.pg_num, self.pg_delay_beg, self.pg_delay_end, self.pg_doppl_beg, self.pg_doppl_end;

    '''
    get the area of channel estimation
    '''
    def getAreaCE(self):
        return self.ce_num, self.ce_delay_beg, self.ce_delay_end, self.ce_doppl_beg, self.ce_doppl_end;

    '''
    check pulse type
    '''
    def isPulseIdeal(self):
        return self.pulse_type == self.PULSE_IDEAL;
    def isPulseRecta(self):
        return self.pulse_type == self.PULSE_RECTA;
    
    '''
    get the pilots matrix
    '''
    def getPilotsMat(self):
        Xp = self.zeros(self.nTimeslotNum, self.nSubcarNum);
        if self.batch_size == self.BATCH_SIZE_NO:
            Xp[self.pk1:self.pk1+self.pk_len, self.pl1:self.pl1+self.pl_len] = self.reshape(self.pilots, self.pk_len, self.pl_len);
        else:
            pilots = np.tile(np.reshape(self.pilots, (self.pk_len, self.pl_len)), (self.batch_size, 1, 1));
            Xp[..., self.pk1:self.pk1+self.pk_len, self.pl1:self.pl1+self.pl_len] = pilots;
        return Xp;
    
    '''
    set content
    @content: a 2D matrix containing pilots, guards and data (if used)
    '''
    def setContent(self, content):
        self.content = content;

    '''
    get content size
    '''
    def getContentSize(self):
        return self.nSubcarNum, self.nTimeslotNum;

    '''
    get content
    @isVector: if true, the returned result is a vector
    '''
    def getContent(self, *, isVector=False):
        # return
        if isVector:
            content = self.reshape(self.content, self.nSubcarNum*self.nTimeslotNum);
        else:
            content = self.content;
        return content;

    '''
    get content - CE (return a matrix)
    '''
    def getContentCE(self):
        if self.ce_num == 0:
            return [];
        else:
            if self.batch_size == self.BATCH_SIZE_NO:
                return self.content[self.ce_doppl_beg:self.ce_doppl_end, self.ce_delay_beg:self.ce_delay_end];
            else:
                return self.content[..., self.ce_doppl_beg:self.ce_doppl_end, self.ce_delay_beg:self.ce_delay_end];
        
    '''
    get content - except CE area (return a vector)
    '''
    def getContentNoCE(self):
        # get data
        if self.ce_num == 0:
            data = self.content.copy();
            data = self.reshape(data, self.nSubcarNum*self.nTimeslotNum);
        else:
            data_num = self.nSubcarNum*self.nTimeslotNum - self.ce_num;
            data = self.zeros(data_num).astype(self.content.dtype);
            data_id = 0;
            for doppl_id in range(self.nTimeslotNum):
                for delay_id in range(self.nSubcarNum):
                    if not self.isInAreaCE(doppl_id, delay_id):
                        if self.batch_size == self.BATCH_SIZE_NO:
                            data[data_id] = self.content[doppl_id, delay_id];
                        else:
                            data[..., data_id] = self.content[..., doppl_id, delay_id];
                        data_id = data_id + 1;
        return data;

    '''
    get content - zero PG area (return a vector)
    '''
    def getContentZeroPG(self):
        # get data
        data = self.content.copy();
        # zero the CE area
        if self.ce_num > 0:
            for doppl_id in range(self.nTimeslotNum):
                for delay_id in range(self.nSubcarNum):
                    if self.isInAreaPG(doppl_id, delay_id):
                        if self.batch_size == self.BATCH_SIZE_NO:
                            data[doppl_id, delay_id] = 0;
                        else:
                            data[..., doppl_id, delay_id] = 0;
        return data;

    '''
    get content - zero CE area (return a vector)
    '''
    def getContentZeroCE(self):
        # get data
        data = self.content.copy();
        # zero the CE area
        if self.ce_num > 0:
            for doppl_id in range(self.nTimeslotNum):
                for delay_id in range(self.nSubcarNum):
                    if self.isInAreaCE(doppl_id, delay_id):
                        if self.batch_size == self.BATCH_SIZE_NO:
                            data[doppl_id, delay_id] = 0;
                        else:
                            data[..., doppl_id, delay_id] = 0;
        return data;
    
    '''
    get content - Data locations
    '''
    def getContentDataLocsMat(self):
        XdLocs = self.ones(self.nTimeslotNum, self.nSubcarNum).astype(bool);
        if self.pilot_type == self.PILOT_TYPE_EM:
            if self.pg_num > 0:
                XdLocs[..., self.pg_doppl_beg:self.pg_doppl_end+1, self.pg_delay_beg:self.pg_delay_end+1] = False;
        return XdLocs;
        
    ###########################################################################
    # private methods
    '''
    validate pilot settings
    '''
    def validP(self):
        # 1st pilot coordinates & pilot length
        if self.pk_len > 0 and self.pl_len > 0:
            if self.pk1 < 0 or self.pk1 >= self.nTimeslotNum:
                raise Exception("Pilot 1st location on the delay axis overflows.");
            if self.pl1 < 0 or self.pl1 >= self.nSubcarNum:
                raise Exception("Pilot 1st location on the Doppler axis overflows.");
            if self.pk1 + self.pk_len - 1 >= self.nTimeslotNum:
                raise Exception("Pilot length on the delay axis overflows.");
            if self.pl1 + self.pl_len - 1 >= self.nSubcarNum:
                raise Exception("Pilot length on the Doppler axis overflows.");
    
    '''
    calculate PG & CE area
    '''
    def calcAreaPGCE(self):
        if self.pl_len > 0 and self.pk_len > 0:
            # overflow check
            if self.pl1 - self.gl_len_neg < 0:
                raise Exception("The guard (neg) on delay axis overflows.");
            if (self.pl1+self.pl_len-1) + self.gl_len_pos >= self.nSubcarNum:
                raise Exception("The guard (pos) on delay axis overflows.");
            if self.pk1 - self.gk_len_neg < 0:
                raise Exception("The guard (neg) on Doppler axis overflows.");
            if (self.pk1+self.pk_len-1) + self.gk_len_pos >= self.nTimeslotNum:
                raise Exception("The guard (pos) on Doppler axis overflows.");
            # calculate PG area
            if self.pilot_type == self.PILOT_TYPE_EM or self.pilot_type == self.PILOT_TYPE_SP:
                # PG area only exist when using embedded pilots
                self.pg_num = ((self.pl_len+self.gl_len_neg+self.gl_len_pos)*(self.pk_len+self.gk_len_neg+self.gk_len_pos)).astype(int);
                self.pg_delay_beg = (self.pl1 - self.gl_len_neg).astype(int);
                self.pg_delay_end = (self.pl1 + self.pl_len - 1 + self.gl_len_pos).astype(int);
                self.pg_doppl_beg = (self.pk1 - self.gk_len_neg).astype(int);
                self.pg_doppl_end = (self.pk1 + self.pk_len - 1 + self.gk_len_pos).astype(int);
            # calulate channel estimate area
            if self.pilot_type == self.PILOT_TYPE_EM or self.pilot_type == self.PILOT_TYPE_SP:
                if self.gl_len_ful:
                    self.ce_delay_beg = 0;
                    self.ce_delay_end = (self.nSubcarNum - 1).astype(int);
                else:
                    self.ce_delay_beg = (self.pg_delay_beg + self.gl_len_neg).astype(int);
                    self.ce_delay_end = (self.pg_delay_end).astype(int);
                if self.gk_len_ful:
                    self.ce_doppl_beg = 0;
                    self.ce_doppl_end = (self.nTimeslotNum - 1).astype(int);
                else:
                    self.ce_doppl_beg = (self.pg_doppl_beg + floor(self.gk_len_neg/2)).astype(int);
                    self.ce_doppl_end = (self.pg_doppl_end - floor(self.gk_len_pos/2)).astype(int);
                self.ce_num = ((self.ce_delay_end - self.ce_delay_beg + 1)*(self.ce_doppl_end - self.ce_doppl_beg + 1)).astype(int);

    '''
    insert data
    @symbols: symbols to map a vector
    '''
    def insertDA(self, symbols):
        symbols = self.squeeze(symbols).copy();
        # input check
        if self.pilot_type == self.PILOT_TYPE_SP:
            data_num = self.nTimeslotNum*self.nSubcarNum - self.zp_len*self.nTimeslotNum;
        else:
            data_num = self.nTimeslotNum*self.nSubcarNum - self.pg_num - self.zp_len*self.nTimeslotNum;
        # input check - symbols
        if not self.isvector(symbols):
            raise Exception("The transmission symbol must be a vector");
        elif symbols.shape[-1] != data_num:
            raise Exception("The transmission symbol number must be %d"%data_num);
        
        # insert
        if data_num == self.nSubcarNum*self.nTimeslotNum:
            self.content = self.reshape(symbols, self.nTimeslotNum, self.nSubcarNum);
        elif data_num == self.nSubcarNum*self.nTimeslotNum - self.zp_len*self.nTimeslotNum:
            self.content = np.append(self.reshape(symbols, self.nTimeslotNum, (self.nSubcarNum - self.zp_len)), self.zeros(self.zp_len, self.nTimeslotNum));
        else:
            symbols_id = 0;
            for doppl_id in range(self.nTimeslotNum):
                for delay_id in range(self.nSubcarNum):
                    if self.isInAreaDA(doppl_id, delay_id):
                        if self.batch_size == self.BATCH_SIZE_NO:
                            self.content[doppl_id, delay_id] = symbols[symbols_id];
                        else:
                            self.content[..., doppl_id, delay_id] = symbols[..., symbols_id];
                        symbols_id = symbols_id + 1;
            assert(symbols_id == data_num);

    '''
    insert pilots
    @pilots_pow: pilot power
    '''
    def insertP(self, pilots_pow):
        # input check
        self.p_len = self.pilots.shape[-1];
        # input check - pilots
        if self.p_len > 0:
            if self.p_len != self.pl_len*self.pk_len:
                raise Exception("The manual pilot input do not have the required number.");
            elif self.p_len > self.nTimeslotNum*self.nSubcarNum:
                raise Exception("The manual pilot input overflows (over the OTFS frame size).");
        # input check - pilot power
        if self.p_len == 0 and pilots_pow == None and self.pk_len>0 and self.pl_len>0:
            raise Exception("The pilots (linear) power is required while no manual pilot input.");
        # initiate pilots if empty
        if self.p_len == 0:
            self.p_len = self.pl_len*self.pk_len;
            if self.p_len > 0:
                self.pilots = sqrt(pilots_pow/2)*(1+1j)*np.ones(self.p_len);
        # allocate pilots
        if self.p_len != 0:
            # allocate pilots
            if self.batch_size == self.BATCH_SIZE_NO:
                self.content[self.pk1:self.pk1+self.pk_len, self.pl1:self.pl1+self.pl_len] = self.content[self.pk1:self.pk1+self.pk_len, self.pl1:self.pl1+self.pl_len] + self.reshape(self.pilots, self.pk_len, self.pl_len);
            else:
                pilots = np.tile(np.reshape(self.pilots, (self.pk_len, self.pl_len)), (self.batch_size, 1, 1));
                self.content[..., self.pk1:self.pk1+self.pk_len, self.pl1:self.pl1+self.pl_len] = self.content[..., self.pk1:self.pk1+self.pk_len, self.pl1:self.pl1+self.pl_len] + pilots;

    '''
    decide whether the current location is in CE area
    @tag:           10->PG, 11->ZP, 12->Data, 20->CE
    @pos_doppl:     the position on the Doppler axis.
    @pos_delay:     the delay position. Not given means the position is for a Doppler-delay vector
    '''
    def isInArea(self, tag, pos_doppl, *args):
        # input check
        if len(args) == 0:
            if pos_doppl < 0 or pos_doppl >= self.nSubcarNum*self.nTimeslotNum:
                raise Exception("The vector position is out of the OTFS size.");
        else:
            if pos_doppl < 0 or pos_doppl >= self.nTimeslotNum:
                raise Exception("The Doppler position is out of the timeslot number.");
            if args[0] < 0 or args[0] > self.nSubcarNum:
                raise Exception("The delay position is out of the subcarrier number.");
        
        # recalculate the position
        if len(args) == 0:
            pos_id = pos_doppl;
            pos_doppl = pos_doppl/self.nSubcarNum;
            if pos_doppl == floor(pos_doppl):
                pos_doppl = floor(pos_doppl);
            else:
                pos_doppl = floor(pos_doppl) + 1;
            delay_pos = pos_id - (pos_doppl-1)*self.nSubcarNum - 1;
        else:
            delay_pos = args[0];
        # decide
        if tag == 10:
            is_in = delay_pos >= self.pg_delay_beg and delay_pos <= self.pg_delay_end and pos_doppl >= self.pg_doppl_beg and pos_doppl <= self.pg_doppl_end;
        elif tag == 11:
            is_in = self.zp_len > 0 and delay_pos >= self.nSubcarNum - self.zp_len;
        elif tag == 12:
            if self.pg_num == 0:
                is_in = True;
            else:
                is_in = delay_pos < self.pg_delay_beg or delay_pos > self.pg_delay_end or pos_doppl < self.pg_doppl_beg or pos_doppl > self.pg_doppl_end;
            if self.zp_len > 0: 
                is_in = is_in and delay_pos < self.nSubcarNum - self.zp_len;
        elif tag == 20:
            is_in = delay_pos >= self.ce_delay_beg and delay_pos <= self.ce_delay_end and pos_doppl >= self.ce_doppl_beg and pos_doppl <= self.ce_doppl_end;
        return is_in;
    
    '''
    estimated the channel using a sngle pilot
    @threshold: the threshold to detect a path
    '''
    def estimateChannel4SingPilot(self, threshold):
        # input check
        if self.p_len != 1:
            raise Exception("There should be only 1 pilot.");
        
        # estimate the channe
        his = [];
        kis = [];
        lis = [];
        for delay_id in range(self.ce_delay_beg, self.ce_delay_end+1):
            for doppl_id in range(self.ce_doppl_beg, self.ce_doppl_end+1):
                pss_ys = np.expand_dims(self.content[doppl_id, delay_id], axis=0) if self.batch_size == self.BATCH_SIZE_NO else self.content[..., doppl_id, delay_id];
                pss_ys_ids_yes = abs(pss_ys) > threshold;
                pss_ys_ids_not = abs(pss_ys) <= threshold;
                li = delay_id - self.pl1;
                ki = doppl_id - self.pk1;
                if self.pulse_type == self.PULSE_IDEAL:
                    pss_beta = exp(-2j*pi*li*ki/self.nSubcarNum/self.nTimeslotNum);
                elif self.pulse_type == self.PULSE_RECTA:
                    pss_beta = exp(2j*pi*self.pl1*ki/self.nSubcarNum/self.nTimeslotNum);
                hi = pss_ys/self.pilots[0]/pss_beta;
                # at least we find one path
                if np.sum(pss_ys_ids_yes, axis=None) > 0:
                    if self.batch_size == self.BATCH_SIZE_NO:
                        his = np.append(his, hi);
                        lis = np.append(lis, li);
                        kis = np.append(kis, ki);
                    else:
                        hi[pss_ys_ids_not] = 0;
                        hi = hi[..., np.newaxis];
                        li = np.tile(li, (self.batch_size, 1));
                        ki = np.tile(ki, (self.batch_size, 1));
                        if isinstance(his, list):
                            his = hi;
                            lis = li;
                            kis = ki;
                        else:
                            his = np.append(his, hi, axis=-1);
                            lis = np.append(lis, li, axis=-1);
                            kis = np.append(kis, ki, axis=-1);
        return his, lis, kis;