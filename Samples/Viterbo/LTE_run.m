close all;
clear;
clc;

%% Load Settings
Settings_LTE_General;
Settings_LTE_Min_UE500;     % UE 500km/h
%Settings_LTE_12x2_UE500;     % UE 500km/h 12x2
%Settings_LTE_12x3_UE500;     % UE 500km/h 12x3
% Settings_LTE_12x4_UE500;     % UE 500km/h 12x4
% Settings_LTE_12x5_UE500;     % UE 500km/h 12x5
% Settings_LTE_12x6_UE500;     % UE 500km/h 12x6
% Settings_LTE_14x12_UE150;   % UE 150km/h
%Settings_LTE_2x2_kmax1;      % UE toy case
Settings_LTE_Files;

%% Show Running Information
c = clock;
fix(c);
fprintf("Current Time %d/%d/%d %d:%d:%.1f\n", c);
fprintf("-------------------------------------------------------------------------\n");
fprintf("|\n");
fprintf("| OTFS Simulation LTE UE %dkm/h\n", UEspeed);
fprintf("| N(time slots) = %d, M(subcarriers) = %d\n", N, M);
fprintf("|\n");
fprintf("-------------------------------------------------------------------------\n");

%% Fixed seed (only for debugging)
% random input bits generation
data_info_bit = randi([0,1],N_bits_perfram,1);
data_temp = bi2de(reshape(data_info_bit,N_syms_perfram,M_bits));
x_origin = qammod(data_temp,M_mod,'gray', 'UnitAveragePower', true);
x = reshape(x_origin,N,M);                  % the input is N*M

% OTFS modulation
s = OTFS_modulation(N,M,x);

% OTFS channel generation
[delay_taps,Doppler_taps,chan_coef] = OTFS_SetChannel(p, lmax, kmax);

% OTFS channel output
[r, noise] = OTFS_PassChannel(N,M,p,delay_taps,Doppler_taps,chan_coef,sigma_2,s, 0);

% OTFS demodulation
y = OTFS_demodulation(N,M,r);

% Generate equivalent channel, Tx and Rx
[~, ~, Heff] = OTFS_DDChannel_Gen(M, N,p,delay_taps,Doppler_taps,chan_coef);
x_DD = x.';
x_DD = x_DD(:);
y_DD = y.';
y_DD = y_DD(:);

% Generate Real values
Heff_real = [real(Heff), -imag(Heff); imag(Heff), real(Heff)];
y_DD_real = [real(y_DD); imag(y_DD)];
x_DD_real = [real(x_DD); imag(x_DD)];