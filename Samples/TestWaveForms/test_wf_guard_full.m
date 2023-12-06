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
N = 9;                          % time slot number
M = 11;                         % subcarrier number
% Gen information symbols (as a column vector)
x_origin_DD = sqrt(1/2)*(ones(N, M) + 1j*ones(N, M));
x_origin_DD(:, 6-2:6+2) = 0;
x_origin_DD(5, 6) = sqrt(1/2)*(1 + 1j);
x_origin = x_origin_DD.';
x_origin = x_origin(:);

% init OTFS
otfs = OTFS(M, N);
% modulate
otfs.modulate(x_origin_DD);
% set the channel
otfs.addChannelPath(1, 1, 1);
H_DD = otfs.getChannel();
% pass the channel
otfs.passChannel(No);
% demodulate
yDD = otfs.demodulate();
Y_DD = reshape(yDD, M, N).';

% calculate the residual
residual = sum(yDD - H_DD*x_origin, "all");
fprintf("The residual is %.16f\n", abs(residual));

% plot
figure("Name", "Full Guards")
subplot(1,2,1)
bar3(abs(x_origin_DD));
title("Before the channel")
subplot(1,2,2)
bar3(abs(Y_DD));
title("After the channel")


