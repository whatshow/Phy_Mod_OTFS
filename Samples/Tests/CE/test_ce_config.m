clear;
clc;
%% setting
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
sympool = qammod(0: M_mod-1, M_mod, 'UnitAveragePower',true);
p = 2;
lmax = 1;
kmax = 1;
guard_delay_num_neg = lmax;
guard_delay_num_pos = lmax;
guard_doppl_num_neg = kmax*2;
guard_doppl_num_pos = kmax*2;
% OTFS & pilots settings
Ns = [8 8 7 7];
Ms = [8 8 7 7];
pilots_num_delays = [2 1 2 1];
pilots_num_doppls = [2 1 2 1];