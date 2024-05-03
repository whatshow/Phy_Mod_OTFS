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
    raise Exception("You have to run matlab script to generate data.");
M = np.squeeze(matlab_data["M"]);
N = np.squeeze(matlab_data["N"]);
pil_pow = np.squeeze(matlab_data["pil_pow"]);
pl_len = np.squeeze(matlab_data["pl_len"]);
pk_len = np.squeeze(matlab_data["pk_len"]);
guard_delay_num_neg = np.squeeze(matlab_data["guard_delay_num_neg"]);
guard_delay_num_pos = np.squeeze(matlab_data["guard_delay_num_pos"]);
xDD = np.squeeze(matlab_data["xDD"]);
yDD_mat = np.squeeze(matlab_data["yDD"]);
content_mat = np.squeeze(matlab_data["content"]);

# build rg
rg = OTFSResGrid(M, N);
rg.setPulse2Ideal();
rg.setPilot2Center(pl_len, pk_len);
rg.setGuard(guard_delay_num_neg, guard_delay_num_pos, guard_doppl_full=True);
rg.map(xDD, pilots_pow=pil_pow);
yDD, his_est, lis_est, kis_est = rg.demap(threshold=1e-10);
content = rg.getContent();

assert(abs(his_est) == 1);
assert(lis_est == [0])
assert(kis_est == [0])

yDD_diff = abs(yDD - yDD_mat);
assert(np.sum(yDD_diff) == 0);
content_diff = abs(content - content_mat);
assert(np.sum(content_diff, axis=None) == 0);

