% 
% test get HDD when using channel estimation (fractional Dopplers)
%
clear;
clc;

%% general configuration
SNR_p = 30; % dB
SNR_d = 10; % dB
No = 0;
pil_pow = 10^(SNR_p/10);
pil_thr = 0;
sig_pow = 10^(SNR_d/10);
p = 2;
lmax = 14;
kmax = 3.5;
% QAM configuration
M_mod = 16;
M_bits = log2(M_mod);
sympool = sqrt(sig_pow)*qammod(0: M_mod-1, M_mod, 'UnitAveragePower',true);

%% test 
nSubcarNum = 256/2;
nTimeslotNum = 32/2;
pilots_num_delay = 1;
pilots_num_doppler = 1;
guard_delay_num_neg = lmax;
guard_delay_num_pos = lmax;
guard_doppl_num_neg = floor(kmax)*2;
guard_doppl_num_pos = floor(kmax)*2;
symbols_len = nTimeslotNum*nSubcarNum-(pilots_num_delay+guard_delay_num_neg+guard_delay_num_pos)*(pilots_num_doppler+guard_doppl_num_neg+guard_doppl_num_pos);
fprintf("Symbol number = %d\n", symbols_len);
fprintf("HDD remove %d columns\n", (pilots_num_delay+guard_delay_num_neg+guard_delay_num_pos)*(pilots_num_doppler+guard_doppl_num_neg+guard_doppl_num_pos));
fprintf("HDD remove %d rows\n", (pilots_num_delay+guard_delay_num_pos)*(pilots_num_doppler+guard_doppl_num_neg/2+guard_doppl_num_pos/2));
fprintf("HDD should be %dx%d\n", (nTimeslotNum*nSubcarNum - (pilots_num_delay+guard_delay_num_pos)*(pilots_num_doppler+guard_doppl_num_neg/2+guard_doppl_num_pos/2)), symbols_len)

otfs = OTFS(nSubcarNum, nTimeslotNum);
otfs.insertPilotsAndGuards(pilots_num_delay, pilots_num_doppler, "pilots_pow", pil_pow, "guard_delay_num_neg", guard_delay_num_neg, "guard_delay_num_pos", guard_delay_num_pos, "guard_doppler_num_neg", guard_doppl_num_neg, "guard_doppler_num_pos", guard_doppl_num_pos);
X_DD_PG = otfs.X_DD;
nbits = randi([0 1], symbols_len*M_bits, 1);
symbols = sqrt(sig_pow)*qammod(nbits, M_mod,'InputType','bit','UnitAveragePower',true);
otfs.modulate(symbols);
X_DD = otfs.X_DD;
otfs.setChannel("p", p, "lmax", lmax, "kmax", kmax);
H_DD = otfs.getChannel();
otfs.passChannel(No);
yDD = otfs.demodulate();
Y_DD = otfs.getYDD();
yDD_diff = abs(yDD - H_DD*symbols);
[gains2_est, delays2_est, dopplers2_est] = otfs.estimateChannel("threshold", 1e-10);
[gains2, delays2, dopplers2] = otfs.getCSI("sort_by_delay_doppler", true);

%% print
fprintf("- the yDD calculation difference is %e\n", sum(yDD_diff)/symbols_len);