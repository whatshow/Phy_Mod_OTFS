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
% OTFS frame
N_syms_perfram = N*M;
N_bits_perfram = N_syms_perfram*M_bits;


%% simulation
% build rg
data_temp = randi([0 1], N_bits_perfram, 1);
xDD = qammod(data_temp,M_mod,'InputType','bit','UnitAveragePower', true);
rg = OTFSResGrid(M, N);
rg.map(xDD);
rg.setPulse2Recta();
%rg.setPulse2Ideal();
% through the channel
otfs = OTFS();
otfs.modulate(rg);
otfs.setChannel(p, lmax, kmax);
otfs.passChannel(No);
[his, lis, kis] = otfs.getCSI("sort_by_delay_doppler", true);
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
