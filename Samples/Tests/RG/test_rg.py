import numpy as np
import sys
import scipy.io
import os
sys.path.append("..");
from OTFSResGrid import OTFSResGrid

# file
project_name = "phy_mod_otfs";
path_folder = os.path.abspath(os.path.dirname(__file__)).lower();
path_folder = path_folder[:path_folder.find(project_name)+len(project_name)];
path_file = os.path.normpath(path_folder+"/_data/Samples/Tests/RG/test_rg.mat");

# load matlab data
try:
    matlab_data = scipy.io.loadmat(path_file);
except FileNotFoundError:
    raise Exception("You have to run `Samples/Tests/test_case_01_gen_data.m` to create data.");
M = np.squeeze(matlab_data["M"]);
N = np.squeeze(matlab_data["N"]);
pil_pow = np.squeeze(matlab_data["pil_pow"]);
pilots_num_delay = np.squeeze(matlab_data["pilots_num_delay"]);
pilots_num_doppl = np.squeeze(matlab_data["pilots_num_doppl"]);
guard_delay_num_neg = np.squeeze(matlab_data["guard_delay_num_neg"]);
guard_delay_num_pos = np.squeeze(matlab_data["guard_delay_num_pos"]);
xDD = np.squeeze(matlab_data["xDD"]);

# build rg
rg = OTFSResGrid(M, N);
rg.setPulse2Ideal();
rg.setPilot2Center(pilots_num_delay, pilots_num_doppl);
rg.setGuard(guard_delay_num_neg, guard_delay_num_pos, guard_doppl_full=True);
rg.map(xDD, pilots_pow=pil_pow);
yDD, his_est, lis_est, kis_est = rg.demap("threshold", 1e-10);