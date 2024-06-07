import numpy as np
import sys
import scipy.io
import os
sys.path.append("..");
from OTFSResGrid import OTFSResGrid
from OTFS import OTFS
from OTFSDetector import OTFSDetector
eps = np.finfo(float).eps;
# batch_size
batch_size = 1;
# file
project_name = "phy_mod_otfs";
path_folder = os.path.abspath(os.path.dirname(__file__)).lower();
path_folder = path_folder[:path_folder.find(project_name)+len(project_name)];
path_file = os.path.normpath(path_folder+"/_data/Samples/Tests/Detector/test_detect_mp_base_full.mat");
# load matlab data
try:
    matlab_data = scipy.io.loadmat(path_file);
except FileNotFoundError:
    raise Exception("You have to run matlab script to generate data.");
M = np.squeeze(matlab_data["M"]);
N = np.squeeze(matlab_data["N"]);
xDD = np.squeeze(matlab_data["xDD"]);
xDD = np.tile(xDD, (batch_size, 1));
yDD_mat = np.squeeze(matlab_data["yDD"]);
yDD_mat = np.tile(yDD_mat, (batch_size, 1));
Y_DD_mat = np.squeeze(matlab_data["Y_DD"]);
Y_DD_mat = np.tile(Y_DD_mat, (batch_size, 1, 1));
# CE
pil_thr = np.squeeze(matlab_data["pil_thr"]);
his = np.squeeze(matlab_data["his"]);
lis = np.squeeze(matlab_data["lis"]);
kis = np.squeeze(matlab_data["kis"]);
his = np.tile(his, (batch_size, 1));
lis = np.tile(lis, (batch_size, 1));
kis = np.tile(kis, (batch_size, 1));
No = np.squeeze(matlab_data["No"]);
noise = np.squeeze(matlab_data["noise"]);
noise = np.tile(noise, (batch_size, 1));
# detect
sympool = np.squeeze(matlab_data["sympool"]);
xDD_est_mat = np.squeeze(matlab_data["xDD_est2"]);
xDD_est_mat = np.tile(xDD_est_mat, (batch_size, 1));

# simulation
# build rg
rg = OTFSResGrid(M, N, batch_size=batch_size);
rg.map(xDD);
rg.setPulse2Recta();
# through the channel
otfs = OTFS(batch_size=batch_size);
otfs.modulate(rg);
otfs.setChannel(his, lis, kis);
otfs.passChannel(noise);
his, lis, kis = otfs.getCSI(sort_by_delay_doppler=True);
rg_rx = otfs.demodulate();
# demapping (CE)
yDD, his_est, lis_est, kis_est = rg_rx.demap(threshold=pil_thr);
Y_DD = rg_rx.getContent();
# detect
od = OTFSDetector(sympool, batch_size=batch_size);
od.useMPBase();
xDD_est = od.detect(rg_rx, his, lis, kis, No);

xDD_est_diff = abs(xDD_est - xDD_est_mat);
assert(np.sum(xDD_est_diff, axis=None)/M/N/batch_size < 1e-13);