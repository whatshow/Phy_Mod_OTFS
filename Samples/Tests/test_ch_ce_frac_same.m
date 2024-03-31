% 
% test get HDD when using channel estimation (fractional Dopplers)
%
clear;
clc;
fprintf("This example shows how to cancel the interference from pilots when fractional Doppler happens.\n");

%% general configuration
SNR_p = 30; % dB
SNR_d = 10; % dB
No = 0;
pil_pow = 10^(SNR_p/10);
pil_thr = 0;
sig_pow = 10^(SNR_d/10);
p = 2;
lmax = 1;
kmax = 1.5;
% QAM configuration
M_mod = 16;
M_bits = log2(M_mod);
sympool = sqrt(sig_pow)*qammod(0: M_mod-1, M_mod, 'UnitAveragePower',true);

%% test 
nSubcarNum = 5;
nTimeslotNum = 5;
pilots_num_delay = 1;
pilots_num_doppler = 1;
guard_delay_num_neg = lmax;
guard_delay_num_pos = lmax;
guard_doppl_num_neg = 2;
guard_doppl_num_pos = 2;
symbols_len = nTimeslotNum*nSubcarNum-(pilots_num_delay+guard_delay_num_neg+guard_delay_num_pos)*(pilots_num_doppler+guard_doppl_num_neg+guard_doppl_num_pos);
symbols_len = nTimeslotNum*nSubcarNum;
fprintf("Symbol number = %d\n", symbols_len);
fprintf("HDD remove %d columns\n", (pilots_num_delay+guard_delay_num_neg+guard_delay_num_pos)*(pilots_num_doppler+guard_doppl_num_neg+guard_doppl_num_pos));
fprintf("HDD remove %d rows\n", (pilots_num_delay+guard_delay_num_pos)*(pilots_num_doppler+guard_doppl_num_neg/2+guard_doppl_num_pos/2));
fprintf("HDD should be %dx%d\n", (nTimeslotNum*nSubcarNum - (pilots_num_delay+guard_delay_num_pos)*(pilots_num_doppler+guard_doppl_num_neg/2+guard_doppl_num_pos/2)), symbols_len)

otfs = OTFS(nSubcarNum, nTimeslotNum);
%otfs.insertPilotsAndGuards(pilots_num_delay, pilots_num_doppler, "pilots_pow", pil_pow, "guard_delay_num_neg", guard_delay_num_neg, "guard_delay_num_pos", guard_delay_num_pos, "guard_doppler_num_neg", guard_doppl_num_neg, "guard_doppler_num_pos", guard_doppl_num_pos);
nbits = randi([0 1], symbols_len*M_bits, 1);
symbols = sqrt(sig_pow)*qammod(nbits, M_mod,'InputType','bit','UnitAveragePower',true);
%symbols = zeros(symbols_len, 1);
otfs.modulate(symbols);
X_DD = otfs.X_DD;
xDD_all = X_DD.';
xDD_all = xDD_all(:); 
otfs.setChannel("p", p, "lmax", lmax, "kmax", kmax);
H_DD = otfs.getChannel();
H_DD_All = otfs.getChannel("only_for_data", false);
otfs.passChannel(No);
yDD = otfs.demodulate();
Y_DD_All = otfs.getYDD();
yDD_all = Y_DD_All.';
yDD_all = yDD_all(:);
yDD_diff = abs(yDD - H_DD*symbols);
yDD_diff_all = abs(yDD_all - H_DD_All*xDD_all);
Y_DD_diff_all = transpose(reshape(yDD_diff_all, nSubcarNum, nTimeslotNum));
%[gains2_est, delays2_est, dopplers2_est] = otfs.estimateChannel("threshold", 1e-10);
% [gains2, delays2, dopplers2] = otfs.getCSI("sort_by_delay_doppler", true);

%% print
fprintf("- the yDD calculation difference is %e\n", sum(yDD_diff)/symbols_len);
fprintf("  - yDD is different because it is influenced by the pilot.\n");
figure("Name", "Y_DD output")
bar3(Y_DD_All);
figure("Name", "Difference in Y_DD")
bar3(Y_DD_diff_all);