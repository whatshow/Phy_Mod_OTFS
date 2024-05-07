import numpy as np
import sys
import scipy.io
import os
sys.path.append("..");
from OTFSResGrid import OTFSResGrid

# settings
pil_pow = 1000;
pil_thr = 1e-10;
N = 8;
M = 8;
lmax = 1;
kmax = 1;
pl_len = 1;
pk_len = 1;
gdn_len = lmax;
gdp_len = lmax;
gkn_len = 2*kmax;
gkp_len = 2*kmax;
# data
xDD_embed = [0.707106781186548 - 0.707106781186548j, \
-0.707106781186548 - 0.707106781186548j, \
0.707106781186548 + 0.707106781186548j, \
-0.707106781186548 - 0.707106781186548j, \
0.707106781186548 - 0.707106781186548j, \
-0.707106781186548 - 0.707106781186548j, \
0.707106781186548 + 0.707106781186548j, \
0.707106781186548 + 0.707106781186548j, \
-0.707106781186548 - 0.707106781186548j, \
0.707106781186548 - 0.707106781186548j, \
0.707106781186548 + 0.707106781186548j, \
0.707106781186548 - 0.707106781186548j, \
0.707106781186548 - 0.707106781186548j, \
0.707106781186548 + 0.707106781186548j, \
0.707106781186548 + 0.707106781186548j, \
0.707106781186548 + 0.707106781186548j, \
-0.707106781186548 + 0.707106781186548j, \
-0.707106781186548 - 0.707106781186548j, \
0.707106781186548 + 0.707106781186548j, \
0.707106781186548 + 0.707106781186548j, \
-0.707106781186548 + 0.707106781186548j, \
0.707106781186548 - 0.707106781186548j, \
-0.707106781186548 + 0.707106781186548j, \
-0.707106781186548 - 0.707106781186548j, \
0.707106781186548 - 0.707106781186548j, \
-0.707106781186548 - 0.707106781186548j, \
0.707106781186548 + 0.707106781186548j, \
-0.707106781186548 + 0.707106781186548j, \
0.707106781186548 + 0.707106781186548j, \
0.707106781186548 + 0.707106781186548j, \
0.707106781186548 + 0.707106781186548j, \
0.707106781186548 - 0.707106781186548j, \
0.707106781186548 - 0.707106781186548j, \
0.707106781186548 + 0.707106781186548j, \
-0.707106781186548 + 0.707106781186548j, \
0.707106781186548 + 0.707106781186548j, \
0.707106781186548 + 0.707106781186548j, \
0.707106781186548 + 0.707106781186548j, \
-0.707106781186548 + 0.707106781186548j, \
0.707106781186548 + 0.707106781186548j];\
xDD_full = [0.707106781186548 - 0.707106781186548j, \
-0.707106781186548 - 0.707106781186548j, \
0.707106781186548 + 0.707106781186548j, \
-0.707106781186548 - 0.707106781186548j, \
0.707106781186548 - 0.707106781186548j, \
-0.707106781186548 - 0.707106781186548j, \
0.707106781186548 + 0.707106781186548j, \
0.707106781186548 + 0.707106781186548j, \
-0.707106781186548 - 0.707106781186548j, \
0.707106781186548 - 0.707106781186548j, \
0.707106781186548 + 0.707106781186548j, \
0.707106781186548 - 0.707106781186548j, \
0.707106781186548 - 0.707106781186548j, \
0.707106781186548 + 0.707106781186548j, \
0.707106781186548 + 0.707106781186548j, \
0.707106781186548 + 0.707106781186548j, \
-0.707106781186548 + 0.707106781186548j, \
-0.707106781186548 - 0.707106781186548j, \
0.707106781186548 + 0.707106781186548j, \
0.707106781186548 + 0.707106781186548j, \
-0.707106781186548 + 0.707106781186548j, \
0.707106781186548 - 0.707106781186548j, \
-0.707106781186548 + 0.707106781186548j, \
-0.707106781186548 - 0.707106781186548j, \
0.707106781186548 - 0.707106781186548j, \
-0.707106781186548 - 0.707106781186548j, \
0.707106781186548 + 0.707106781186548j, \
-0.707106781186548 - 0.707106781186548j, \
0.707106781186548 - 0.707106781186548j, \
-0.707106781186548 - 0.707106781186548j, \
0.707106781186548 + 0.707106781186548j, \
0.707106781186548 + 0.707106781186548j, \
-0.707106781186548 - 0.707106781186548j, \
0.707106781186548 - 0.707106781186548j, \
0.707106781186548 + 0.707106781186548j, \
0.707106781186548 - 0.707106781186548j, \
0.707106781186548 - 0.707106781186548j, \
0.707106781186548 + 0.707106781186548j, \
0.707106781186548 + 0.707106781186548j, \
0.707106781186548 + 0.707106781186548j, \
-0.707106781186548 + 0.707106781186548j, \
-0.707106781186548 - 0.707106781186548j, \
0.707106781186548 + 0.707106781186548j, \
0.707106781186548 + 0.707106781186548j, \
-0.707106781186548 + 0.707106781186548j, \
0.707106781186548 - 0.707106781186548j, \
-0.707106781186548 + 0.707106781186548j, \
-0.707106781186548 - 0.707106781186548j, \
0.707106781186548 - 0.707106781186548j, \
-0.707106781186548 - 0.707106781186548j, \
0.707106781186548 + 0.707106781186548j, \
-0.707106781186548 + 0.707106781186548j, \
0.707106781186548 + 0.707106781186548j, \
0.707106781186548 + 0.707106781186548j, \
0.707106781186548 + 0.707106781186548j, \
0.707106781186548 - 0.707106781186548j, \
0.707106781186548 - 0.707106781186548j, \
0.707106781186548 + 0.707106781186548j, \
-0.707106781186548 + 0.707106781186548j, \
0.707106781186548 + 0.707106781186548j, \
0.707106781186548 + 0.707106781186548j, \
0.707106781186548 + 0.707106781186548j, \
-0.707106781186548 + 0.707106781186548j, \
0.707106781186548 + 0.707106781186548j]; \
batch_size = 10;

'''
rg - embedded, ideal, full guard
'''
# create data
# batch 1
xDD_embed_batch_1 = np.tile(xDD_embed, (1, 1));
xDD_embed_batch_1_row = np.expand_dims(xDD_embed_batch_1, axis=-2);
xDD_embed_batch_1_col = np.expand_dims(xDD_embed_batch_1, axis=-1);
xDD_embed_batch_1_list = xDD_embed_batch_1.tolist();
xDD_embed_batch_1_row_list = xDD_embed_batch_1_row.tolist();
xDD_embed_batch_1_col_list = xDD_embed_batch_1_col.tolist();
# batch n
xDD_embed_batch_n = np.tile(xDD_embed, (batch_size, 1));
for b_id in range(batch_size):
    xDD_embed_batch_n[b_id, ...] = np.roll(xDD_embed_batch_n[b_id, ...], b_id, axis=-1); 
xDD_embed_batch_n_row = np.expand_dims(xDD_embed_batch_n, axis=-2);
xDD_embed_batch_n_col = np.expand_dims(xDD_embed_batch_n, axis=-1);
xDD_embed_batch_n_list = xDD_embed_batch_n.tolist();
xDD_embed_batch_n_row_list = xDD_embed_batch_n_row.tolist();
xDD_embed_batch_n_col_list = xDD_embed_batch_n_col.tolist();
# X_DD
X_DD_embed = np.reshape(xDD_embed, (8,5));
X_DD_embed_batch_1 = np.reshape(xDD_embed_batch_1, (1, 8,5));
X_DD_embed_batch_n = np.reshape(xDD_embed_batch_n, (batch_size, 8,5));
# to inputs and outputs
embed_xDDs = [xDD_embed, 
             xDD_embed_batch_1, 
             xDD_embed_batch_1_row, 
             xDD_embed_batch_1_col,
             xDD_embed_batch_1_list, 
             xDD_embed_batch_1_row_list, 
             xDD_embed_batch_1_col_list,
             xDD_embed_batch_n,
             xDD_embed_batch_n_row,
             xDD_embed_batch_n_col,
             xDD_embed_batch_n_list,
             xDD_embed_batch_n_row_list,
             xDD_embed_batch_n_col_list];
embed_X_DDs = [X_DD_embed,
               X_DD_embed_batch_1,
               X_DD_embed_batch_1,
               X_DD_embed_batch_1,
               X_DD_embed_batch_1,
               X_DD_embed_batch_1,
               X_DD_embed_batch_1,
               X_DD_embed_batch_n,
               X_DD_embed_batch_n,
               X_DD_embed_batch_n,
               X_DD_embed_batch_n,
               X_DD_embed_batch_n,
               X_DD_embed_batch_n];
embed_batch = [None,
               1,1,1,1,1,1,
               batch_size, batch_size, batch_size, batch_size, batch_size, batch_size];
# test
for case_id in range(len(embed_batch)):
    if embed_batch[case_id] is None:
        rg = OTFSResGrid(M, N);
    else:
        rg = OTFSResGrid(M, N, batch_size=embed_batch[case_id]);
    rg.setPulse2Ideal();
    rg.setPilot2Center(pl_len, pk_len);
    rg.setGuard(gdn_len, gdp_len, guard_doppl_full=True);
    rg.map(embed_xDDs[case_id], pilots_pow=pil_pow);
    yDD, his_est, lis_est, kis_est = rg.demap(threshold=1e-10);
    Y_DD,_,_,_ = rg.demap(isDataVec=False, threshold=1e-10);
    Y_DD = np.delete(Y_DD, [2,3,4], axis=-1);
    content = rg.getContent();
    if embed_batch[case_id] is None:
        contentD = np.reshape(yDD, (8,6));
    else:
        contentD = np.reshape(yDD, (embed_batch[case_id], 8,6));
    contentD = np.delete(contentD, 2, axis=-1);
    if embed_batch[case_id] is None:
        assert(abs(his_est) == 1);
    else:
        assert(abs(np.sum(his_est)) == embed_batch[case_id]);
    assert(np.sum(lis_est) == 0);
    assert(np.sum(kis_est) ==0);
    assert(abs(np.sum(contentD - embed_X_DDs[case_id])) == 0);
    assert(abs(np.sum(Y_DD - embed_X_DDs[case_id])) == 0);

'''
rg - superimposed, ideal, full guard
'''
# create data
# batch 1
xDD_full_batch_1 = np.tile(xDD_full, (1, 1));
xDD_full_batch_1_row = np.expand_dims(xDD_full_batch_1, axis=-2);
xDD_full_batch_1_col = np.expand_dims(xDD_full_batch_1, axis=-1);
xDD_full_batch_1_list = xDD_full_batch_1.tolist();
xDD_full_batch_1_row_list = xDD_full_batch_1_row.tolist();
xDD_full_batch_1_col_list = xDD_full_batch_1_col.tolist();
# batch n
xDD_full_batch_n = np.tile(xDD_full, (batch_size, 1));
for b_id in range(batch_size):
    xDD_full_batch_n[b_id, ...] = np.roll(xDD_full_batch_n[b_id, ...], b_id, axis=-1); 
xDD_full_batch_n_row = np.expand_dims(xDD_full_batch_n, axis=-2);
xDD_full_batch_n_col = np.expand_dims(xDD_full_batch_n, axis=-1);
xDD_full_batch_n_list = xDD_full_batch_n.tolist();
xDD_full_batch_n_row_list = xDD_full_batch_n_row.tolist();
xDD_full_batch_n_col_list = xDD_full_batch_n_col.tolist();
# X_DD
X_DD_full = np.reshape(xDD_full, (8,8));
X_DD_full_batch_1 = np.reshape(xDD_full_batch_1, (1, 8,8));
X_DD_full_batch_n = np.reshape(xDD_full_batch_n, (batch_size, 8,8));
# to inputs and outputs
full_xDDs = [xDD_full, 
             xDD_full_batch_1, 
             xDD_full_batch_1_row, 
             xDD_full_batch_1_col,
             xDD_full_batch_1_list, 
             xDD_full_batch_1_row_list, 
             xDD_full_batch_1_col_list,
             xDD_full_batch_n,
             xDD_full_batch_n_row,
             xDD_full_batch_n_col,
             xDD_full_batch_n_list,
             xDD_full_batch_n_row_list,
             xDD_full_batch_n_col_list];
full_X_DDs = [X_DD_full,
               X_DD_full_batch_1,
               X_DD_full_batch_1,
               X_DD_full_batch_1,
               X_DD_full_batch_1,
               X_DD_full_batch_1,
               X_DD_full_batch_1,
               X_DD_full_batch_n,
               X_DD_full_batch_n,
               X_DD_full_batch_n,
               X_DD_full_batch_n,
               X_DD_full_batch_n,
               X_DD_full_batch_n];
full_batch = [None,
               1,1,1,1,1,1,
               batch_size, batch_size, batch_size, batch_size, batch_size, batch_size];
# test
for case_id in range(len(full_batch)):
    if full_batch[case_id] is None:
        rg = OTFSResGrid(M, N);
    else:
        rg = OTFSResGrid(M, N, batch_size=full_batch[case_id]);
    rg.setPulse2Ideal();
    rg.setPilot2SuperImposed();
    rg.setPilot2Center(pl_len, pk_len);
    rg.setGuard(gdn_len, gdp_len, guard_doppl_full=True);
    rg.map(full_xDDs[case_id], pilots_pow=pil_pow);
    yDD, his_est, lis_est, kis_est = rg.demap(threshold=1e-10);
    Y_DD,_,_,_ = rg.demap(isDataVec=False, threshold=1e-10);
    if full_batch[case_id] is None:
        Y_DD_rebuild = np.reshape(yDD, (8, 8));
    else:
        Y_DD_rebuild = np.reshape(yDD, (full_batch[case_id], 8, 8));
    Y_DD_pilot = (1+1j)*np.sqrt(pil_pow/2);
    if full_batch[case_id] is not None:
        Y_DD_pilot = Y_DD_pilot*full_batch[case_id];
    assert(abs(np.sum(Y_DD_rebuild - full_X_DDs[case_id]) - Y_DD_pilot) < 1e-13);
    assert(abs(np.sum(Y_DD - full_X_DDs[case_id]) - Y_DD_pilot) < 1e-13);