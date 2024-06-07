import numpy as np
import sys
import scipy.io
import os
sys.path.append("../../..");
# from OTFSResGrid import OTFSResGrid
# from OTFS import OTFS
from whatshow_phy_mod_otfs import OTFS, OTFSResGrid, OTFSDetector
import matplotlib.pyplot as plt

batch_size = 16;

# # file
# project_name = "phy_mod_otfs";
# path_folder = os.path.abspath(os.path.dirname(__file__)).lower();
# path_folder = path_folder[:path_folder.find(project_name)+len(project_name)];
# path_file = os.path.normpath(path_folder+"/_data/Samples/Tests/RG/test_rg.mat");

# # load matlab data
# try:
#     matlab_data = scipy.io.loadmat(path_file);
# except FileNotFoundError:
#     raise Exception("You have to run matlab script to generate data.");
# M = np.squeeze(matlab_data["M"]);
# N = np.squeeze(matlab_data["N"]);
# pil_pow = np.squeeze(matlab_data["pil_pow"]);
# pl_len = np.squeeze(matlab_data["pl_len"]);
# pk_len = np.squeeze(matlab_data["pk_len"]);
# guard_delay_num_neg = np.squeeze(matlab_data["guard_delay_num_neg"]);
# guard_delay_num_pos = np.squeeze(matlab_data["guard_delay_num_pos"]);
# xDD = np.squeeze(matlab_data["xDD"]);
# xDD = np.tile(xDD, (batch_size, 1));
# yDD_mat = np.squeeze(matlab_data["yDD"]);
# yDD_mat = np.tile(yDD_mat, (batch_size, 1));
# content_mat = np.squeeze(matlab_data["content"]);
# content_mat = np.tile(content_mat, (batch_size, 1, 1));
 
M=N=8
batch_size=128
pl_len = pk_len =1
guard_delay_num_neg=2
guard_delay_num_pos=2


# SNR_d = 10
SNR_p = 20
arr_snr_d = np.array([10,12,14,16,18,20])
mse = []

for SNR_d in arr_snr_d:
    sigma2 = 1/np.power(10,SNR_d/10)
     
    pil_pow = sigma2 *np.power(10,SNR_p/10)
    
    constel = [-0.7071-0.7071j, -0.7071+0.7071j, 0.7071-0.7071j, 0.7071+0.7071j];
    
    # no batch
    # X_DD
    xDD_idx = np.random.randint(4, size=(batch_size, M*N-40));
    xDD = np.take(constel,xDD_idx);
    
    # build rg
    rg = OTFSResGrid(M, N, batch_size=batch_size);
    rg.setPulse2Recta();
    rg.setPilot2Center(pl_len, pk_len);
    rg.setGuard(guard_delay_num_neg, guard_delay_num_pos, guard_doppl_full=True);
    
    # rg.setPilot2SuperImposed()
    rg.map(xDD, pilots_pow=pil_pow);
    
    otfs = OTFS(batch_size=batch_size);
    otfs.modulate(rg);
    p = 6;
    lmax = 2;
    kmax = 1;
     
    
    otfs.setChannel(p, lmax, kmax, force_frac=True);
    otfs.passChannel(sigma2);
    his, lis, kis = otfs.getCSI();
    HDD = otfs.getChannel(data_only=False);
    
    rg_rx = otfs.demodulate();
    
     
    yDD, his_est, lis_est, kis_est = rg_rx.demap(threshold=3*np.sqrt(sigma2));
    
    # _, his_est2, lis_est2, kis_est2 = rg_rx.demap(isData=False, threshold=1e-10);
    
    HDD_est = otfs.getChannel(his=his_est, lis=lis_est, kis=kis_est, data_only=False);
    MSE_Heff = np.mean(np.square(abs(HDD - HDD_est))); 
    mse.append(MSE_Heff)
    Y_DD = rg_rx.getContent();
    yDD = rg_rx.getContent(isVector=True);
    
    
plt.plot(arr_snr_d,mse)
plt.semilogy()
plt.show()



# yDD_diff = abs(yDD - yDD_mat);
# assert(np.sum(yDD_diff) == 0);
# content_diff = abs(content - content_mat);
# assert(np.sum(content_diff, axis=None) == 0);

