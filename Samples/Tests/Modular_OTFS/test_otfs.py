import numpy as np
import sys
import scipy.io
import os
sys.path.append("../../..");
from OTFS import OTFS
from OTFSResGrid import OTFSResGrid

# other settings
M = 8;
N = 8;
pl_len = 1;
pk_len = 1;
guard_delay_num_neg = 2;
guard_delay_num_pos = 2;
pil_pow = 1;
p = 2;
lmax = 2;
kmax = 3;
# data
constel = [-0.7071-0.7071j, -0.7071+0.7071j, 0.7071-0.7071j, 0.7071+0.7071j];
# data - iterative
case_num = 3;
batch_sizes = [None, 1, 5];

# getChannel
# getChannel - embedded
for i in range(case_num):
    batch_size = batch_sizes[i];
    # Tx
    if batch_size is None:
        xDD = np.take(constel, np.random.randint(4, size=(M*N-40)));
    else:
        xDD = np.take(constel, np.random.randint(4, size=(batch_size, M*N-40)));
        xDD = xDD[..., np.newaxis];
    if batch_size is None:
        rg = OTFSResGrid(M, N);
    else:
        rg = OTFSResGrid(M, N, batch_size=batch_size);
    rg.setPulse2Ideal();
    rg.setPilot2Center(pl_len, pk_len);
    rg.setGuard(guard_delay_num_neg, guard_delay_num_pos, guard_doppl_full=True);
    rg.map(xDD, pilots_pow=pil_pow);
    # channel
    if batch_size is None:
        otfs = OTFS();
    else:
        otfs = OTFS(batch_size=batch_size);
    otfs.modulate(rg);
    otfs.setChannel(p, lmax, kmax, force_frac=True);
    otfs.passChannel(0);
    his, lis, kis = otfs.getCSI();
    HDD = otfs.getChannel();
    rg_rx = otfs.demodulate();
    yDD, his_est, lis_est, kis_est = rg_rx.demap(threshold=1e-10);
    HDD_est = otfs.getChannel(his_est, lis_est, kis_est);
    HDD_diff = abs(HDD - HDD_est);
    if batch_size is None:
        yDD_diff = abs(yDD - HDD@xDD);
    else:
        yDD_diff = abs(yDD - np.squeeze(HDD@xDD, axis=-1));
    if batch_size is None:
        assert(np.sum(HDD_diff)/M/N < 1e-13);
        assert(np.sum(yDD_diff)/M/N < 1e-13);
    else:
        assert(np.sum(HDD_diff)/M/N/batch_size < 1e-13);
        assert(np.sum(yDD_diff)/M/N/batch_size < 1e-13);
# getChannel - superimposed
for i in range(case_num):
    batch_size = batch_sizes[i];
    # Tx
    if batch_size is None:
        xDD = np.take(constel, np.random.randint(4, size=(M*N)));
    else:
        xDD = np.take(constel, np.random.randint(4, size=(batch_size, M*N)));
        xDD = xDD[..., np.newaxis];
    if batch_size is None:
        rg = OTFSResGrid(M, N);
    else:
        rg = OTFSResGrid(M, N, batch_size=batch_size);
    rg.setPulse2Ideal();
    rg.setPilot2SuperImposed();                                             # superimposed
    rg.setPilot2Center(pl_len, pk_len);
    rg.setGuard(guard_delay_num_neg, guard_delay_num_pos, guard_doppl_full=True);
    rg.map(xDD, pilots_pow=pil_pow);
    # channel
    if batch_size is None:
        otfs = OTFS();
    else:
        otfs = OTFS(batch_size=batch_size);
    otfs.modulate(rg);
    otfs.setChannel(p, lmax, kmax, force_frac=True);
    otfs.passChannel(0);
    his, lis, kis = otfs.getCSI();
    HDD = otfs.getChannel(data_only=False);
    rg_rx = otfs.demodulate();
    _, his_est, lis_est, kis_est = rg_rx.demap(threshold=1e-10);
    HDD_est = otfs.getChannel(his_est, lis_est, kis_est, data_only=False);
    
    # reget xDD and yDD
    if batch_size is None:
        xDD_all = rg.getContent(isVector=True);
    else:
        xDD_all = rg.getContent(isVector=True);
        xDD_all = xDD_all[..., np.newaxis];
    yDD_all = rg_rx.getContent(isVector=True);
    if batch_size is None:
        yDD_diff = abs(yDD_all - HDD@xDD_all);
    else:
        yDD_diff = abs(yDD_all - np.squeeze(HDD@xDD_all, axis=-1));
    if batch_size is None:
        assert(np.sum(yDD_diff)/M/N < 1e-13);
    else:
        assert(np.sum(yDD_diff)/M/N/batch_size < 1e-13);
    # channel estimation check
    Y_DD = rg_rx.getContent();
    if batch_size is None:
        for pid in range(len(his)):
            hi = his_est[pid];
            li = int(lis_est[pid]);
            ki = int(kis_est[pid]);
                
            residual = Y_DD[rg.pk1+ki, rg.pl1+li] - hi*(1+1j)*np.sqrt(pil_pow/2)*np.exp(-2j*np.pi*li*ki/M/N);
            assert(abs(residual) < 1e-13);
    else:
        for batch_id in range(batch_size):
            for pid in range(len(his[batch_id,:])):
                hi = his_est[batch_id, pid];
                li = int(lis_est[batch_id, pid]);
                ki = int(kis_est[batch_id, pid]);
                residual = Y_DD[batch_id, rg.pk1+ki, rg.pl1+li] - hi*(1+1j)*np.sqrt(pil_pow/2)*np.exp(-2j*np.pi*li*ki/M/N);
                assert(abs(residual) < 1e-13);
    