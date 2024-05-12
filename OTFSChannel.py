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
    # private methods
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
        