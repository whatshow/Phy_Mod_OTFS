clear;
clc;

%% settings
SNR_p = 30; % dB
SNR_d = 18; % dB
No = 10^(-SNR_p/10);
pil_pow = 10^(SNR_p/10);
sig_pow = 10^(SNR_d/10);
pil_thr = 3*sqrt(No);
% OTFS configuration
N = 2;                          % time slot number
M = 2;                          % subcarrier number
% channel
p = 2;
lmax = 1;
kmax = 1;
% QAM configuration
M_mod = 4;
M_bits = log2(M_mod);
sympool = qammod(0: M_mod-1, M_mod, 'UnitAveragePower',true);

%% simulation
% build rg
xDD = [0.707106781186548 + 0.707106781186548i;...
-0.707106781186548 + 0.707106781186548i;...
0.707106781186548 + 0.707106781186548i;...
0.707106781186548 + 0.707106781186548i];
rg = OTFSResGrid(M, N);
rg.map(xDD);
rg.setPulse2Recta();
%rg.setPulse2Ideal();
% through the channel
otfs = OTFS();
otfs.modulate(rg);
his = [0.5,-0.384327912144509+0.210675546668676i];
lis = [0,1];
kis = [0,-1];
otfs.setChannel(his, lis, kis);
noise = [-0.0368960842935724 + 0.0114770685416059i;...
-0.0134216174854600 + 0.0360584963075490i;...
0.00277687537275079 - 0.0105179801136721i;...
-0.0343638173157366 - 0.00324561290652126i];
otfs.passChannel(noise);
rg_rx = otfs.demodulate();
% demapping (CE)
[yDD, ~, ~, ~] = rg_rx.demap("threshold", pil_thr);
Y_DD = rg_rx.getContent();
% detect
od = OTFSDetector(sympool);
od.useMPBase();
xDD_est = od.detect(rg_rx, his, lis, kis, No);
X_DD_est2 = OTFS_mp_detector(N,M,M_mod,p,lis,kis,his,No,Y_DD, "constellation", sympool).';
xDD_est2 = X_DD_est2(:);
assert(sum(xDD_est2-xDD_est) == 0);


%% File
path_folder = "./_data/Samples/Tests/Detector/";
path_file = path_folder + "test_detect_mp_base_full.mat";
% create the folder if not exist
if ~exist(path_folder, 'dir')
    mkdir(path_folder);
end
% delete the file if exist
if exist(path_file, 'file')
    delete(path_file)
end
save(path_file);
