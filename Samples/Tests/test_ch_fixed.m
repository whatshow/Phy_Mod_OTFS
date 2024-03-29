clear;
clc;

% QAM configuration
M_mod = 4;                                                                  % size of constellation
M_bits = log2(M_mod);
sympool = qammod(0: M_mod-1, M_mod, 'UnitAveragePower',true);               % Generate the symbol pool
sympool_real = unique(real(sympool));
SNR = 10; % dB
%No = 1/10^(SNR/10); % linear
No = 0;

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
otfs.setChannel("delays", [0, 1], "Dopplers", [2, 3], "gains", [0.5, 0.5]);
H_DD = otfs.getChannel();
% pass the channel
r = otfs.passChannel(No);
% demodulate
yDD = otfs.demodulate();

% calculate the residual
residual = sum(yDD - H_DD*x_origin, "all");
fprintf("The residual is %.16f\n", abs(residual));