clear;
clc;

%% File
path_folder = "./_data/Samples/Tests/OTFS/";
path_file = path_folder + "test_otfs_inte.mat";
% create the folder if not exist
if ~exist(path_folder, 'dir')
    mkdir(path_folder);
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
N = 5;
M = 7;
lmax = 1;
kmax = 1;
pl_len = 1;
pk_len = 1;
guard_delay_num_neg = lmax;
guard_delay_num_pos = lmax;
% 
xDD = [-0.707106781186548 + 0.707106781186548i, ...
-0.707106781186548 + 0.707106781186548i, ...
-0.707106781186548 + 0.707106781186548i, ...
-0.707106781186548 - 0.707106781186548i, ...
-0.707106781186548 - 0.707106781186548i, ...
-0.707106781186548 + 0.707106781186548i, ...
-0.707106781186548 + 0.707106781186548i, ...
0.707106781186548 + 0.707106781186548i, ...
0.707106781186548 - 0.707106781186548i, ...
0.707106781186548 - 0.707106781186548i, ...
0.707106781186548 - 0.707106781186548i, ...
0.707106781186548 - 0.707106781186548i, ...
-0.707106781186548 + 0.707106781186548i, ...
0.707106781186548 - 0.707106781186548i, ...
-0.707106781186548 - 0.707106781186548i, ...
-0.707106781186548 + 0.707106781186548i, ...
-0.707106781186548 + 0.707106781186548i, ...
-0.707106781186548 + 0.707106781186548i, ...
-0.707106781186548 - 0.707106781186548i, ...
0.707106781186548 - 0.707106781186548i];
%% build rg
rg1 = OTFSResGrid(M, N);
rg1.setPulse2Ideal();
rg1.setPilot2Center(pl_len, pk_len);
rg1.setGuard(guard_delay_num_neg, guard_delay_num_pos, "guard_doppl_full", true);
rg1.map(xDD, "pilots_pow", pil_pow);
rg2 = OTFSResGrid(M, N);
rg2.setPulse2Recta();
rg2.setPilot2Center(pl_len, pk_len);
rg2.setGuard(guard_delay_num_neg, guard_delay_num_pos, "guard_doppl_full", true);
rg2.map(xDD, "pilots_pow", pil_pow);
%% pass the channel
otfs = OTFS();
otfs.setChannel([0.5,0.1463-0.5938j], [0, 1], [1, -1]);
otfs.modulate(rg1);
otfs.passChannel(0);
rg1_rx = otfs.demodulate();
otfs.modulate(rg2);
otfs.passChannel(0);
rg2_rx = otfs.demodulate();
%% demap
[yDD1, his_est1, lis_est1, kis_est1] = rg1_rx.demap("threshold", 1e-10);
content1 = rg1_rx.getContent();
[yDD2, his_est2, lis_est2, kis_est2] = rg2_rx.demap("threshold", 1e-10);
content2 = rg2_rx.getContent();
%% save demap result
if ~exist(path_file, 'file')
    save(path_file);
end