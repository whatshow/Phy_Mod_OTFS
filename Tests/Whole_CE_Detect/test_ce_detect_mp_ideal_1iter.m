clear;
clc;

%% settings
SNR_p = 30; % dB
SNR_d = 18; % dB
No = 1;
pil_pow = 10^(SNR_p/10);
sig_pow = 10^(SNR_d/10);
pil_thr = 3*sqrt(No);
% OTFS configuration
N = 64;                          % time slot number
M = 16;                          % subcarrier number
% channel
p = 3;
lmax = 1;
kmax = 1;
% pilot settings
pilots_num_delay = 1;
pilots_num_doppl = 1;
guard_delay_num_neg = lmax;
guard_delay_num_pos = lmax;
guard_doppl_num_neg = kmax*2;
guard_doppl_num_pos = kmax*2;
% QAM configuration
M_mod = 16;
M_bits = log2(M_mod);
sympool = sqrt(sig_pow)*qammod(0: M_mod-1, M_mod, 'UnitAveragePower',true);
% OTFS frame
N_syms_perfram = N*M-(pilots_num_delay+guard_delay_num_neg+guard_delay_num_pos)*N;
N_bits_perfram = N_syms_perfram*M_bits;


%% simulation
% build rg
data_temp = randi([0 1], N_bits_perfram, 1);
xDD = sqrt(sig_pow)*qammod(data_temp,M_mod,'InputType','bit','UnitAveragePower', true);
rg = OTFSResGrid(M, N);
rg.setPilot2Center(pilots_num_delay, pilots_num_doppl);
rg.setGuard(guard_delay_num_neg, guard_delay_num_pos, 'guard_doppl_full', true);
rg.map(xDD, "pilots_pow", pil_pow);
%rg.setPulse2Recta();
rg.setPulse2Ideal();
% through the channel
otfs = OTFS();
otfs.modulate(rg);
otfs.setChannel(p, lmax, kmax);
otfs.passChannel(No);
[his, lis, kis] = otfs.getCSI("sort_by_delay_doppler", true);
rg_rx = otfs.demodulate();
% demapping (CE)
[~, his_est, lis_est, kis_est] = rg_rx.demap("isData", false, "threshold", pil_thr);
% detect
od = OTFSDetector(sympool);
od.useMPBase();
xDD_est = od.detect(rg_rx, his_est, lis_est, kis_est, No);
SER = (abs(xDD_est - xDD) > eps)/N_syms_perfram;
fprintf("MP differenct is %e\n", sum(SER));

