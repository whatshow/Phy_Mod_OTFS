import numpy as np
import sys
import scipy.io
import os
sys.path.append("..");
from OTFSResGrid import OTFSResGrid

batch_size = 16;

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
xDD = np.tile(xDD, (batch_size, 1));
yDD_mat = np.squeeze(matlab_data["yDD"]);
yDD_mat = np.tile(yDD_mat, (batch_size, 1));
content_mat = np.squeeze(matlab_data["content"]);
content_mat = np.tile(content_mat, (batch_size, 1, 1));

# build rg
rg = OTFSResGrid(M, N, batch_size=batch_size);
rg.setPulse2Ideal();
rg.setPilot2Center(pl_len, pk_len);
rg.setGuard(guard_delay_num_neg, guard_delay_num_pos, guard_doppl_full=True);
rg.map(xDD, pilots_pow=pil_pow);
yDD, his_est, lis_est, kis_est = rg.demap(threshold=1e-10);
content = rg.getContent();

assert(sum(abs(his_est)) == batch_size);
assert(sum(lis_est) == 0)
assert(sum(kis_est) == 0)

yDD_diff = abs(yDD - yDD_mat);
assert(np.sum(yDD_diff) == 0);
content_diff = abs(content - content_mat);
assert(np.sum(content_diff, axis=None) == 0);

