clear;
clc;

%% settings
SNR_p = 30; % dB
SNR_d = 10; % dB
No = 1;
pil_pow = 10^(SNR_p/10);
pil_thr = 3*sqrt(1/pil_pow);
sig_pow = 10^(SNR_d/10);
% OTFS configuration
N = 16;                          % time slot number
M = 16;                          % subcarrier number
% channel
p = 4;
lmax = 2;
kmax = 2;
% pilot settings
pilots_num_delay = 1;
pilots_num_doppler = 1;
guard_delay_num_neg = lmax;
guard_delay_num_pos = lmax;
guard_doppl_num_neg = kmax*2;
guard_doppl_num_pos = kmax*2;
% QAM configuration
M_mod = 16;
M_bits = log2(M_mod);
sympool = sqrt(sig_pow)*qammod(0: M_mod-1, M_mod, 'UnitAveragePower',true);
% OTFS frame
N_syms_perfram = N*M-(pilots_num_delay+guard_delay_num_neg+guard_delay_num_pos)*(pilots_num_doppler+guard_doppl_num_neg+guard_doppl_num_pos);
N_bits_perfram = N_syms_perfram*M_bits;


%% simulation
% Gen information symbols (as a column vector)
data_temp = randi([0 1], N_bits_perfram, 1);
x_origin = sqrt(sig_pow)*qammod(data_temp,M_mod,'InputType','bit','UnitAveragePower', true);

% init OTFS
otfs = OTFS(M, N, "pilot_loc_type", OTFS.PILOT_LOC_CENTER);
% insert 
otfs.insertPilotsAndGuards(pilots_num_delay, pilots_num_doppler, "pilots_pow", pil_pow, "guard_delay_num_neg", guard_delay_num_neg, "guard_delay_num_pos", guard_delay_num_pos, "guard_doppler_num_neg", guard_doppl_num_neg, "guard_doppler_num_pos", guard_doppl_num_pos);
% modulate
otfs.modulate(x_origin);
% set the channel
otfs.setChannel("p", p, "lmax", lmax, "kmax", kmax);
% pass the channel
r = otfs.passChannel(No);
% demodulate
yDD = otfs.demodulate();
YDD = otfs.getYDD();
% estimate the channel
[chan_coef_est, delay_taps_est, doppler_taps_est] = otfs.estimateChannel("threshold", pil_thr);
[chan_coef, delay_taps, doppler_taps] = otfs.getCSI("sort_by_delay_doppler", true);
taps = length(chan_coef_est);
% detect
xDD_est = otfs.detect(OTFS.DETECT_MP_BASE, OTFS.DETECT_CSI_CE, No, sympool, "sym_map", true);
xDD_est_perfect = otfs.detect(OTFS.DETECT_MP_BASE, OTFS.DETECT_CSI_PERFECT, No, sympool, "sym_map", true);

% SER cal
SER_Perfect = sum(xDD_est_perfect ~= x_origin)/N_syms_perfram;
fprintf("SER (Perfect) %e\n", SER_Perfect);
SER_CE = sum(xDD_est ~= x_origin)/N_syms_perfram;
fprintf("SER (CE) %e\n", SER_CE);
