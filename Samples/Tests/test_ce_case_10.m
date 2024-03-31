% 
% test illegal pilot allocation
%
clear;
clc;

% OTFS settings
nSubcarNum = 256/2;
nTimeslotNum = 32/2;
% channel
p = 6;
kmax = 3.5;
lmax = 14;
% pilot & data
SNR_d = 10;                                 % data SNR range (dB)
SNR_p = 10;                                 % pilot SNR (dB) (we have the condition of)
No = 1;                                     % noise power (linear)
% pilot
pil_pow = 10.^(SNR_p./10);
pilots_num_delay = 1;
pilots_num_doppler = 1;
guard_delay_num_neg = lmax;
guard_delay_num_pos = lmax;
guard_doppl_num_neg = ceil(kmax)*2;
guard_doppl_num_pos = ceil(kmax)*2;
% signal
sig_pow = 10.^(SNR_d./10);
sig_len = nSubcarNum*nTimeslotNum - (pilots_num_delay+guard_delay_num_neg+guard_delay_num_pos)*(pilots_num_doppler+guard_doppl_num_neg+guard_doppl_num_pos);
% Tx Settings
M_mod = 16;
M_bits = log2(M_mod);
sympool_norm = qammod(0: M_mod-1, M_mod, 'UnitAveragePower',true);
% init otfs
otfs = OTFS(nSubcarNum, nTimeslotNum, "pilot_loc_type", OTFS.PILOT_LOC_CENTER);
% insert pilots
disp("The guard allocation on the Doppler axis is going to overflow.");
otfs.insertPilotsAndGuards(pilots_num_delay, pilots_num_doppler, "pilots_pow", pil_pow, "guard_delay_num_neg", guard_delay_num_neg, "guard_delay_num_pos", guard_delay_num_pos, "guard_doppler_num_neg", guard_doppl_num_neg, "guard_doppler_num_pos", guard_doppl_num_pos);
% xDD_syms
nbits = randi([0,1],sig_len*M_bits,1);
xDD_syms = sqrt(sig_pow/2)*qammod(nbits, M_mod,'InputType','bit','UnitAveragePower',true);
% modulate
otfs.modulate(xDD_syms);