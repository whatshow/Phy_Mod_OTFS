% 
% test all pilots allocation for **center pilots and reduced guard** (no noise)
%
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

%% test center - reduced guard - even size - even pilot number
N = 8;
M = 8;
pilots_num_delay = 2;
pilots_num_doppler = 2;
guard_delay_num_neg = lmax;
guard_delay_num_pos = lmax;
guard_doppl_num_neg = kmax*2;
guard_doppl_num_pos = kmax*2;
symbols_len = N*M-(pilots_num_delay+guard_delay_num_neg+guard_delay_num_pos)*(pilots_num_doppler+guard_doppl_num_neg+guard_doppl_num_pos);
otfs = OTFS(M, N, "pilot_type", OTFS.PILOT_SINGLE_SISO, "pilot_loc_type", OTFS.PILOT_LOC_CENTER, "GUARD_TYPE", OTFS.GUARD_REDUCED);
otfs.insertPilotsAndGuards(pilots_num_delay, pilots_num_doppler, "pilots_pow", pil_pow, "guard_delay_num_neg", guard_delay_num_neg, "guard_delay_num_pos", guard_delay_num_pos, "guard_doppler_num_neg", guard_doppl_num_neg, "guard_doppler_num_pos", guard_doppl_num_pos);
X_DD_PG1 = otfs.X_DD;
nbits = randi([0 1], symbols_len*M_bits, 1);
symbols = sqrt(sig_pow)*qammod(nbits, M_mod,'InputType','bit','UnitAveragePower',true);
otfs.modulate(symbols);
X_DD1 = otfs.X_DD;
otfs.setChannel("p", p, "lmax", lmax, "kmax", kmax);
H_DD1 = otfs.getChannel();
otfs.passChannel(No);

%% test center - reduced guard - even size - odd pilot number
N = 8;
M = 8;
pilots_num_delay = 1;
pilots_num_doppler = 1;
guard_delay_num_neg = lmax;
guard_delay_num_pos = lmax;
guard_doppl_num_neg = kmax*2;
guard_doppl_num_pos = kmax*2;
symbols_len = N*M-(pilots_num_delay+guard_delay_num_neg+guard_delay_num_pos)*(pilots_num_doppler+guard_doppl_num_neg+guard_doppl_num_pos);
otfs = OTFS(M, N, "pilot_type", OTFS.PILOT_SINGLE_SISO, "pilot_loc_type", OTFS.PILOT_LOC_CENTER, "GUARD_TYPE", OTFS.GUARD_REDUCED);
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
[gains2_est, delays2_est, dopplers2_est] = otfs.estimateChannel("threshold", 1e-10);
[gains2, delays2, dopplers2] = otfs.getCSI();

%% test center - reduced guard - odd size - even pilot number
N = 7;
M = 7;
pilots_num_delay = 2;
pilots_num_doppler = 2;
guard_delay_num_neg = lmax;
guard_delay_num_pos = lmax;
guard_doppl_num_neg = kmax*2;
guard_doppl_num_pos = kmax*2;
symbols_len = N*M-(pilots_num_delay+guard_delay_num_neg+guard_delay_num_pos)*(pilots_num_doppler+guard_doppl_num_neg+guard_doppl_num_pos);
otfs = OTFS(M, N, "pilot_type", OTFS.PILOT_SINGLE_SISO, "pilot_loc_type", OTFS.PILOT_LOC_CENTER, "GUARD_TYPE", OTFS.GUARD_REDUCED);
otfs.insertPilotsAndGuards(pilots_num_delay, pilots_num_doppler, "pilots_pow", pil_pow, "guard_delay_num_neg", guard_delay_num_neg, "guard_delay_num_pos", guard_delay_num_pos, "guard_doppler_num_neg", guard_doppl_num_neg, "guard_doppler_num_pos", guard_doppl_num_pos);
X_DD_PG3 = otfs.X_DD;
nbits = randi([0 1], symbols_len*M_bits, 1);
symbols = sqrt(sig_pow)*qammod(nbits, M_mod,'InputType','bit','UnitAveragePower',true);
otfs.modulate(symbols);
X_DD3 = otfs.X_DD;
otfs.setChannel("p", p, "lmax", lmax, "kmax", kmax);
H_DD3 = otfs.getChannel();
otfs.passChannel(No);

%% test center - reduced guard - odd size - odd pilot number
N = 7;
M = 7;
pilots_num_delay = 1;
pilots_num_doppler = 1;
guard_delay_num_neg = lmax;
guard_delay_num_pos = lmax;
guard_doppl_num_neg = kmax*2;
guard_doppl_num_pos = kmax*2;
symbols_len = N*M-(pilots_num_delay+guard_delay_num_neg+guard_delay_num_pos)*(pilots_num_doppler+guard_doppl_num_neg+guard_doppl_num_pos);
otfs = OTFS(M, N, "pilot_type", OTFS.PILOT_SINGLE_SISO, "pilot_loc_type", OTFS.PILOT_LOC_CENTER, "GUARD_TYPE", OTFS.GUARD_REDUCED);
otfs.insertPilotsAndGuards(pilots_num_delay, pilots_num_doppler, "pilots_pow", pil_pow, "guard_delay_num_neg", guard_delay_num_neg, "guard_delay_num_pos", guard_delay_num_pos, "guard_doppler_num_neg", guard_doppl_num_neg, "guard_doppler_num_pos", guard_doppl_num_pos);
X_DD_PG4 = otfs.X_DD;
nbits = randi([0 1], symbols_len*M_bits, 1);
symbols = sqrt(sig_pow)*qammod(nbits, M_mod,'InputType','bit','UnitAveragePower',true);
otfs.modulate(symbols);
X_DD4 = otfs.X_DD;
otfs.setChannel("p", p, "lmax", lmax, "kmax", kmax);
H_DD4 = otfs.getChannel();
otfs.passChannel(No);

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
figure("name", "pilots & guards & datda");
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