% 
% test all pilots allocation for **center pilots and reduced guard** (no noise)
%
clear;
clc;

%% general configuration
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
sympool = sqrt(sig_pow)*qammod(0: M_mod-1, M_mod, 'UnitAveragePower',true);
p = 2;
lmax = 1;
kmax = 1;

%% test CE
fprintf("test CE\n");

%% test CE - reduced guard
fprintf("  - reduced guard\n");
%% test CE - reduced guard - even size, even pilot number
fprintf("    - even size, even pilot number\n");
N = 8;
M = 8;
pilots_num_delay = 2;
pilots_num_doppler = 2;
guard_delay_num_neg = lmax;
guard_delay_num_pos = lmax;
guard_doppl_num_neg = kmax*2;
guard_doppl_num_pos = kmax*2;
symbols_len = N*M-(pilots_num_delay+guard_delay_num_neg+guard_delay_num_pos)*(pilots_num_doppler+guard_doppl_num_neg+guard_doppl_num_pos);
otfs = OTFS(M, N, "pilot_loc_type", OTFS.PILOT_LOC_CENTER);
otfs.insertPilotsAndGuards(pilots_num_delay, pilots_num_doppler, "pilots_pow", pil_pow, "guard_delay_num_neg", guard_delay_num_neg, "guard_delay_num_pos", guard_delay_num_pos, "guard_doppler_num_neg", guard_doppl_num_neg, "guard_doppler_num_pos", guard_doppl_num_pos);
X_DD_PG1 = otfs.X_DD;
nbits = randi([0 1], symbols_len*M_bits, 1);
symbols = sqrt(sig_pow)*qammod(nbits, M_mod,'InputType','bit','UnitAveragePower',true);
otfs.modulate(symbols);
X_DD1 = otfs.X_DD;
otfs.setChannel("p", p, "lmax", lmax, "kmax", kmax);
H_DD1 = otfs.getChannel();
otfs.passChannel(No);
otfs.demodulate();
Y_DD1 = otfs.getYDD();
[gains1_est, delays1_est, dopplers1_est] = otfs.estimateChannel("threshold", 1e-10);
[gains1, delays1, dopplers1] = otfs.getCSI("sort_by_delay_doppler", true);

%% test CE - reduced guard - even size - odd pilot number
fprintf("    - even size, odd pilot number\n");
N = 8;
M = 8;
pilots_num_delay = 1;
pilots_num_doppler = 1;
guard_delay_num_neg = lmax;
guard_delay_num_pos = lmax;
guard_doppl_num_neg = kmax*2;
guard_doppl_num_pos = kmax*2;
symbols_len = N*M-(pilots_num_delay+guard_delay_num_neg+guard_delay_num_pos)*(pilots_num_doppler+guard_doppl_num_neg+guard_doppl_num_pos);
otfs = OTFS(M, N, "pilot_loc_type", OTFS.PILOT_LOC_CENTER);
otfs.insertPilotsAndGuards(pilots_num_delay, pilots_num_doppler, "pilots_pow", pil_pow, "guard_delay_num_neg", guard_delay_num_neg, "guard_delay_num_pos", guard_delay_num_pos, "guard_doppler_num_neg", guard_doppl_num_neg, "guard_doppler_num_pos", guard_doppl_num_pos);
X_DD_PG2 = otfs.X_DD;
nbits = randi([0 1], symbols_len*M_bits, 1);
symbols = sqrt(sig_pow)*qammod(nbits, M_mod,'InputType','bit','UnitAveragePower',true);
otfs.modulate(symbols);
X_DD2 = otfs.X_DD;
otfs.setChannel("p", p, "lmax", lmax, "kmax", kmax);
H_DD2 = otfs.getChannel();
otfs.passChannel(No);
otfs.demodulate();
Y_DD2 = otfs.getYDD();
[gains2_est, delays2_est, dopplers2_est] = otfs.estimateChannel("threshold", 1e-10);
[gains2, delays2, dopplers2] = otfs.getCSI("sort_by_delay_doppler", true);
fprintf("      - the CE difference\n");
fprintf("        - gain: %e, delay: %e, Doppler: %e\n", sum(abs(gains2 - gains2_est)), sum(delays2-delays2_est), sum(dopplers2_est - dopplers2));


%% test CE - reduced guard - odd size - even pilot number
fprintf("    - odd size, even pilot number\n");
N = 7;
M = 7;
pilots_num_delay = 2;
pilots_num_doppler = 2;
guard_delay_num_neg = lmax;
guard_delay_num_pos = lmax;
guard_doppl_num_neg = kmax*2;
guard_doppl_num_pos = kmax*2;
symbols_len = N*M-(pilots_num_delay+guard_delay_num_neg+guard_delay_num_pos)*(pilots_num_doppler+guard_doppl_num_neg+guard_doppl_num_pos);
otfs = OTFS(M, N, "pilot_loc_type", OTFS.PILOT_LOC_CENTER);
otfs.insertPilotsAndGuards(pilots_num_delay, pilots_num_doppler, "pilots_pow", pil_pow, "guard_delay_num_neg", guard_delay_num_neg, "guard_delay_num_pos", guard_delay_num_pos, "guard_doppler_num_neg", guard_doppl_num_neg, "guard_doppler_num_pos", guard_doppl_num_pos);
X_DD_PG3 = otfs.X_DD;
nbits = randi([0 1], symbols_len*M_bits, 1);
symbols = sqrt(sig_pow)*qammod(nbits, M_mod,'InputType','bit','UnitAveragePower',true);
otfs.modulate(symbols);
X_DD3 = otfs.X_DD;
otfs.setChannel("p", p, "lmax", lmax, "kmax", kmax);
H_DD3 = otfs.getChannel();
otfs.passChannel(No);
otfs.demodulate();
Y_DD3 = otfs.getYDD();
[gains3_est, delays3_est, dopplers3_est] = otfs.estimateChannel("threshold", 1e-10);
[gains3, delays3, dopplers3] = otfs.getCSI("sort_by_delay_doppler", true);

%% test center - reduced guard - odd size - odd pilot number
fprintf("    - odd size, odd pilot number\n");
N = 7;
M = 7;
pilots_num_delay = 1;
pilots_num_doppler = 1;
guard_delay_num_neg = lmax;
guard_delay_num_pos = lmax;
guard_doppl_num_neg = kmax*2;
guard_doppl_num_pos = kmax*2;
symbols_len = N*M-(pilots_num_delay+guard_delay_num_neg+guard_delay_num_pos)*(pilots_num_doppler+guard_doppl_num_neg+guard_doppl_num_pos);
otfs = OTFS(M, N, "pilot_loc_type", OTFS.PILOT_LOC_CENTER);
otfs.insertPilotsAndGuards(pilots_num_delay, pilots_num_doppler, "pilots_pow", pil_pow, "guard_delay_num_neg", guard_delay_num_neg, "guard_delay_num_pos", guard_delay_num_pos, "guard_doppler_num_neg", guard_doppl_num_neg, "guard_doppler_num_pos", guard_doppl_num_pos);
X_DD_PG4 = otfs.X_DD;
nbits = randi([0 1], symbols_len*M_bits, 1);
symbols = sqrt(sig_pow)*qammod(nbits, M_mod,'InputType','bit','UnitAveragePower',true);
otfs.modulate(symbols);
X_DD4 = otfs.X_DD;
otfs.setChannel("p", p, "lmax", lmax, "kmax", kmax);
H_DD4 = otfs.getChannel();
otfs.passChannel(No);
otfs.demodulate();
Y_DD4 = otfs.getYDD();
[gains4_est, delays4_est, dopplers4_est] = otfs.estimateChannel("threshold", 1e-10);
[gains4, delays4, dopplers4] = otfs.getCSI("sort_by_delay_doppler", true);
fprintf("      - the CE difference\n");
fprintf("        - gain: %e, delay: %e, Doppler: %e\n", sum(abs(gains4 - gains4_est)), sum(delays4-delays4_est), sum(dopplers4_est - dopplers4));

%% plot
% pilots & guards
figure("name", "pilots & guards");
subplot(2,2,1);
bar3(abs(X_DD_PG1));
title("reduced guard - even size - even pilot number")
subplot(2,2,2);
bar3(abs(X_DD_PG2));
title("reduced guard - even size - odd pilot number")
subplot(2,2,3);
bar3(abs(X_DD_PG3));
title("reduced guard - odd size - even pilot number")
subplot(2,2,4);
bar3(abs(X_DD_PG4));
title("reduced guard - odd size - odd pilot number")
% pilots & guards & data
figure("name", "pilots & guards & data");
subplot(2,2,1);
bar3(abs(X_DD1));
title("reduced guard - even size - even pilot number")
subplot(2,2,2);
bar3(abs(X_DD2));
title("reduced guard - even size - odd pilot number")
subplot(2,2,3);
bar3(abs(X_DD3));
title("reduced guard - odd size - even pilot number")
subplot(2,2,4);
bar3(abs(X_DD4));
title("reduced guard - odd size - odd pilot number")
% pilots & guards & data
figure("name", "After the channel");
subplot(2,2,1);
bar3(abs(Y_DD1));
title("reduced guard - even size - even pilot number")
subplot(2,2,2);
bar3(abs(Y_DD2));
title("reduced guard - even size - odd pilot number")
subplot(2,2,3);
bar3(abs(Y_DD3));
title("reduced guard - odd size - even pilot number")
subplot(2,2,4);
bar3(abs(Y_DD4));
title("reduced guard - odd size - odd pilot number")