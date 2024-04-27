clear;
clc;

%% File
path_folder = "./_data/Samples/Tests/RG/";
path_file = path_folder + "test_rg.mat";
% create the folder if not exist
if ~exist(path_folder, 'dir')
    mkdir(path_folder);
end
% delete the file if exist
if exist(path_file, 'file')
    delete(path_file)
end

%% settings
SNR_p = 30; % dB
SNR_d = 10; % dB
No = 0;
% No = 1;
pil_pow = 10^(SNR_p/10);
pil_thr = 0;
% pil_thr = 3*sqrt(1/pil_pow);
sig_pow = 10^(SNR_d/10);
M_mod = 4;
M_bits = log2(M_mod);
%
N = 8;
M = 8;
lmax = 1;
kmax = 1;
pl_len = 1;
pk_len = 1;
guard_delay_num_neg = lmax;
guard_delay_num_pos = lmax;
guard_doppl_num_neg = kmax*2;
guard_doppl_num_pos = kmax*2;
% 
symbols_len = N*M-(pl_len+guard_delay_num_neg+guard_delay_num_pos)*N;
nbits = randi([0 1], symbols_len*M_bits, 1);
xDD = qammod(nbits, M_mod,'InputType','bit','UnitAveragePower',true);
%% save config to file
save(path_file);
%% build rg
rg = OTFSResGrid(M, N);
rg.setPulse2Ideal();
rg.setPilot2Center(pl_len, pk_len);
rg.setGuard(guard_delay_num_neg, guard_delay_num_pos, "guard_doppl_full", true);
rg.map(xDD, "pilots_pow", pil_pow);
[yDD, his_est, lis_est, kis_est] = rg.demap("threshold", 1e-10);
content = rg.getContent();
%% save demap result
save(path_file, "yDD", "content", "-append");