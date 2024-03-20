% MP detection test using fixed data
%
% This script is to test MP detector performs exact the same as the method 
% from Interference Cancellation and Iterative Detection for Orthogonal 
% Time Frequency Space Modulation by Raviteja Patchava, Yi Hong, and 
% Emanuele Viterbo in 2018.

clear;
clc;

% QAM configuration
M_mod = 4;                                                                  % size of constellation
M_bits = log2(M_mod);
sympool = qammod(0: M_mod-1, M_mod, 'UnitAveragePower',true);               % Generate the symbol pool
SNR = 20; % dB
No = 1/10^(SNR/10); % linear

% OTFS configuration
N = 8;                          % time slot number
M = 8;                         % subcarrier number
N_syms_perfram = N*M;           % number of symbols per frame
N_bits_perfram = N*M*M_bits;    % number of bits per frame

% iterative detect
rng(1)
N_fram = 10^4;
SERs_MP_OTFS_Fram = zeros(N_fram,1);
SERs_MP_Vite_Fram = zeros(N_fram,1);
parfor ifram = 1:N_fram
    % Gen information symbols (as a column vector)
    data_temp = randi([0,1],N_bits_perfram,1);
    x_origin = qammod(data_temp, M_mod,'InputType','bit','UnitAveragePower',true);
    % init OTFS
    otfs = OTFS(M, N);
    % modulate
    otfs.modulate(x_origin);
    % set the channel
    otfs.setChannel("p", 4, "lmax", 3, "kmax", 3);
    % pass the channel
    r = otfs.passChannel(No);
    % demodulate
    yDD = otfs.demodulate();
    % get CSI
    taps = 4;
    [chan_coef,delay_taps,doppler_taps]  = otfs.getCSI();
    % detect
    xDD_est = otfs.detect(OTFS.DETECT_MP_BASE, OTFS.DETECT_CSI_PERFECT, No, sympool);
    X_DD_est_viterbo = OTFS_mp_detector(N,M,M_mod,taps,delay_taps,doppler_taps,chan_coef,No, otfs.getYDD(), "constellation", sympool);
    xDD_est_viterbo = X_DD_est_viterbo.';
    xDD_est_viterbo = xDD_est_viterbo(:);
    
    if sum(xDD_est_viterbo ~= xDD_est) > eps
        error("The proposed MP in OTFS is wrong.");
    end
    
    SERs_MP_OTFS_Fram(ifram) = sum(xDD_est - x_origin > eps)/N/M;
    SERs_MP_Vite_Fram(ifram) = sum(xDD_est_viterbo - x_origin > eps)/N/M;
end

SERs_MP_Vite = mean(SERs_MP_Vite_Fram);
fprintf("Viterbo's MP SER is %.16f\n", SERs_MP_Vite);
SERs_MP_OTFS = mean(SERs_MP_OTFS_Fram);
fprintf("Our MP SER is %.16f\n", SERs_MP_Vite);

% difference
SERs_diff = abs(SERs_MP_OTFS - SERs_MP_Vite);
fprintf("The SER difference between Viterbo's base MP and ours is %.16f\n", SERs_diff);

