import numpy as np
import sys
import scipy.io
import os
sys.path.append("..");
from OTFSResGrid import OTFSResGrid
from OTFS import OTFS
eps = np.finfo(float).eps;

batch_size = 16;

# file
project_name = "phy_mod_otfs";
path_folder = os.path.abspath(os.path.dirname(__file__)).lower();
path_folder = path_folder[:path_folder.find(project_name)+len(project_name)];
path_file = os.path.normpath(path_folder+"/_data/Samples/Tests/OTFS/test_otfs_inte.mat");
# load matlab data
try:
    matlab_data = scipy.io.loadmat(path_file);
except FileNotFoundError:
    raise Exception("You have to run matlab script to generate data.");
M = np.squeeze(matlab_data["M"]);
N = np.squeeze(matlab_data["N"]);
pil_pow = np.squeeze(matlab_data["pil_pow"]);
pl_len = np.squeeze(matlab_data["pl_len"]);
pk_len = np.squeeze(matlab_data["pk_len"]);
guard_delay_num_neg = np.squeeze(matlab_data["guard_delay_num_neg"]);
guard_delay_num_pos = np.squeeze(matlab_data["guard_delay_num_pos"]);
xDD = np.tile(np.squeeze(matlab_data["xDD"]), (batch_size, 1));
yDD1_mat = np.tile(np.squeeze(matlab_data["yDD1"]), (batch_size, 1));
yDD2_mat = np.tile(np.squeeze(matlab_data["yDD2"]), (batch_size, 1));
content1_mat = np.tile(np.squeeze(matlab_data["content1"]), (batch_size, 1, 1));
content2_mat = np.tile(np.squeeze(matlab_data["content2"]), (batch_size, 1, 1));
# CE
his_est1_mat = np.squeeze(matlab_data["his_est1"]);
lis_est1_mat = np.squeeze(matlab_data["lis_est1"]);
kis_est1_mat = np.squeeze(matlab_data["kis_est1"]);
his_est2_mat = np.squeeze(matlab_data["his_est2"]);
lis_est2_mat = np.squeeze(matlab_data["lis_est2"]);
kis_est2_mat = np.squeeze(matlab_data["kis_est2"]);


# build rg
rg1 = OTFSResGrid(M, N, batch_size=batch_size);
rg1.setPulse2Ideal();
rg1.setPilot2Center(pl_len, pk_len);
rg1.setGuard(guard_delay_num_neg, guard_delay_num_pos, guard_doppl_full=True);
rg1.map(xDD, pilots_pow=pil_pow);
rg2 = OTFSResGrid(M, N, batch_size=batch_size);
rg2.setPulse2Recta();
rg2.setPilot2Center(pl_len, pk_len);
rg2.setGuard(guard_delay_num_neg, guard_delay_num_pos, guard_doppl_full=True);
rg2.map(xDD, pilots_pow=pil_pow);
# pass the channel
otfs = OTFS(batch_size=batch_size);
otfs.modulate(rg1);
his = np.tile([0.5,0.1463-0.5938j], (batch_size, 1));
lis = np.tile([0, 1], (batch_size, 1));
kis = np.tile([1, -1], (batch_size, 1));
otfs.setChannel(his, lis, kis);
otfs.passChannel(0);
rg1_rx = otfs.demodulate();
otfs.modulate(rg2);
otfs.passChannel(0);
rg2_rx = otfs.demodulate();
# demap
yDD1, his_est1, lis_est1, kis_est1 = rg1_rx.demap(threshold=1e-10);
content1 = rg1_rx.getContent();
# check
yDD1_diff = abs(yDD1 - yDD1_mat);
assert(np.sum(yDD1_diff)/M/N/batch_size < 1e-13);
content1_diff = abs(content1 - content1_mat);
assert(np.sum(content1_diff, axis=None)/M/N/batch_size < 1e-13);
assert(np.sum(abs(his_est1 - his_est1_mat)/batch_size, axis=None)< 1e-13);

yDD2, his_est2, lis_est2, kis_est2 = rg2_rx.demap(threshold=1e-10);
content2 = rg2_rx.getContent();
# check
yDD2_diff = abs(yDD2 - yDD2_mat);
assert(np.sum(yDD2_diff)/M/N/batch_size < 1e-13);
content2_diff = abs(content2 - content2_mat);
assert(np.sum(content2_diff, axis=None)/M/N/batch_size <1e-13);
assert(np.sum(abs(his_est2 - his_est2_mat)/batch_size, axis=None)< 1e-13);

