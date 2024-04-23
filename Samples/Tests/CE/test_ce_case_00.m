% Case 00:  Overflow on the Doppler axis
%           
% In this example, we need 8+8+1=17 grids on the Doppler axis while we 
% only have 16 grids. Therefore, you should observe an exception.
%
clear;
clc;
%% OTFS settings
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
pilots_num_doppl = 1;
guard_delay_num_neg = lmax;
guard_delay_num_pos = lmax;
guard_doppl_num_neg = 8;
guard_doppl_num_pos = 8;
% signal
sig_pow = 10.^(SNR_d./10);
sig_len = nSubcarNum*nTimeslotNum - (pilots_num_delay+guard_delay_num_neg+guard_delay_num_pos)*(pilots_num_doppl+guard_doppl_num_neg+guard_doppl_num_pos);
% Tx Settings
M_mod = 16;
M_bits = log2(M_mod);
sympool_norm = qammod(0: M_mod-1, M_mod, 'UnitAveragePower',true);

%% Symbols
% xDD_syms
nbits = randi([0,1],sig_len*M_bits,1);
xDD_syms = sqrt(sig_pow/2)*qammod(nbits, M_mod,'InputType','bit','UnitAveragePower',true);

%% RG
rg = OTFSResGrid(nSubcarNum, nTimeslotNum);
rg.setPilot2Center(pilots_num_delay, pilots_num_doppl);
rg.setGuard(guard_delay_num_neg, guard_delay_num_pos, guard_doppl_num_neg, guard_doppl_num_pos);
rg.map(xDD_syms, "pilots_pow", pil_pow);