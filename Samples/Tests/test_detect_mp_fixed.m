% MP detection test using fixed data
%
% This script is to test MP detector performs exact the same as the method 
% from Interference Cancellation and Iterative Detection for Orthogonal 
% Time Frequency Space Modulation by Raviteja Patchava, Yi Hong, and 
% Emanuele Viterbo in 2018.

clear;
clc;

% QAM configuration
M_mod = 16;                                                                  % size of constellation
M_bits = log2(M_mod);
sympool = qammod(0: M_mod-1, M_mod, 'UnitAveragePower',true);               % Generate the symbol pool
SNR = 10; % dB
No = 1/10^(SNR/10); % linear
%No = 0;

% OTFS configuration
N = 4;                          % time slot number
M = 4;                         % subcarrier number
N_syms_perfram = N*M;           % number of symbols per frame
N_bits_perfram = N*M*M_bits;    % number of bits per frame


% Gen information symbols (as a column vector)
data_temp = [3;3;2;3;0;2;2;0;1;2;1;1;3;3;2;0];
x_origin = qammod(data_temp,M_mod,'gray', 'UnitAveragePower', true);

% init OTFS
otfs = OTFS(M, N);
% modulate
otfs.modulate(x_origin);
% set the channel
taps = 2;
chan_coef = [0.5, 0.5];
delay_taps = [0, 1];
doppler_taps = [2, 3];
otfs.setChannel("gains", chan_coef, "delays", delay_taps, "Dopplers", doppler_taps);
% pass the channel
r = otfs.passChannel(No);
% demodulate
yDD = otfs.demodulate();
% detect
xDD_est = otfs.detect(OTFS.DETECT_MP_BASE, OTFS.DETECT_CSI_PERFECT, No, sympool);
X_DD_est_viterbo = OTFS_mp_detector(N,M,M_mod,taps,delay_taps,doppler_taps,chan_coef,No, otfs.getYDD(), "constellation", sympool);
xDD_est_viterbo = X_DD_est_viterbo.';
xDD_est_viterbo = xDD_est_viterbo(:);

% difference
xDD_est_diff = abs(xDD_est - xDD_est_viterbo);
fprintf("The difference between Viterbo's base MP and ours is %.16f\n", sum(xDD_est_diff));

