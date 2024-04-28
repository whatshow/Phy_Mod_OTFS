import numpy as np
import sys
import scipy.io
import os
sys.path.append("..");
from OTFSResGrid import OTFSResGrid
from OTFS import OTFS
from OTFSDetector import OTFSDetector
eps = np.finfo(float).eps;
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
yDD_mat = np.squeeze(matlab_data["yDD"]);
Y_DD_mat = np.squeeze(matlab_data["Y_DD"]);
# CE
pil_thr = np.squeeze(matlab_data["pil_thr"]);
his = np.squeeze(matlab_data["his"]);
lis = np.squeeze(matlab_data["lis"]);
kis = np.squeeze(matlab_data["kis"]);
No = np.squeeze(matlab_data["No"]);
noise = np.squeeze(matlab_data["noise"]);
# detect
sympool = np.squeeze(matlab_data["sympool"]);
xDD_est_mat = np.squeeze(matlab_data["xDD_est2"]);

# simulation
# build rg
rg = OTFSResGrid(M, N);
rg.map(xDD);
rg.setPulse2Recta();
# through the channel
otfs = OTFS();
otfs.modulate(rg);
otfs.setChannel(his, lis, kis);
otfs.passChannel(noise);
his, lis, kis = otfs.getCSI(sort_by_delay_doppler=True);
rg_rx = otfs.demodulate();
# demapping (CE)
yDD, his_est, lis_est, kis_est = rg_rx.demap(threshold=pil_thr);
Y_DD = rg_rx.getContent();
# detect
od = OTFSDetector(sympool);
od.useMPBase();
xDD_est = od.detect(rg_rx, his, lis, kis, No);

xDD_est_diff = abs(xDD_est - xDD_est_mat);
assert(np.sum(xDD_est_diff)/M/N < 1e-13);