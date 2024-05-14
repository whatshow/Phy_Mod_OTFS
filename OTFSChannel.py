import numpy as np
from whatshow_toolbox import MatlabFuncHelper


class OTFSChannel(MatlabFuncHelper):
    # channel coefficients
    config = None;
    his = None;
    lis = None;
    kis = None;
    
    ###########################################################################
    # General Methods
    '''
    constructor
    @config:                the configuration of the entire OTFS sysmtem
    @batch_size(obsolete):  the batch size
    '''
    def __init__(self, *, config=None, batch_size=None):
        if config is not None:
            self.config = config;
        if batch_size is not None:
            self.batch_size = batch_size;
    
    '''
    
    '''
    def CSI2HisMat(self, his, lis, kis):
        pass
    def HisMat2CSI(self, his, lis, kis):
        pass
    
    ###########################################################################
    # obsolete methods
    def obsCSIList2Mat(his, lis, kis, lmax, kmin, kmax):
        pass
    
    def obsCSIMat2List(his_mat, lmax, kmin, kmax):
        his = [];
        lis = [];
        kis = [];
        
        return his, lis, kis;