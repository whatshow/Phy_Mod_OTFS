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
fft_size = 128;                 % fft size
N = 9;                          % time slot number
M = 11;                         % subcarrier number
tap_pos = [1, 1; 2, 1; 1, 11];
descriptions = ["Base"; "1 Doppler"; "11 Delay"];

% generante data
sigs = zeros(fft_size*N, length(descriptions));
sigs_base = zeros(M*N, length(descriptions));
for id = 1:length(descriptions)
    % get configuration
    tap_pos_x = tap_pos(id, 1);
    tap_pos_y = tap_pos(id, 2);
    
    X_DD = zeros(N, M);
    X_DD(tap_pos_x, tap_pos_y) = sqrt(1/2)*(1 + 1j);
    
    % init OTFS
    otfs = OTFS(M, N);
    % modulate
    otfs.modulate(X_DD);
    % get the symbol in the time domain
    sigs(:, id) = otfs.getS("fft_size", fft_size);
    sigs_base(:, id) = otfs.getS("fft_size", M);
end

%% plot (high resolution)
figure("Name", "Waveform (high resolution)")
xindices = 1:fft_size:(N*fft_size+1);
% base
subplot(3,2,1)
plot(real(sigs(:, 1)));
title("Base (real)");
grid on;
xlim([1, N*fft_size]);
xticks(xindices);
subplot(3,2,2)
plot(imag(sigs(:, 1)));
title("Base (imag)");
grid on;
xlim([1, N*fft_size]);
xticks(xindices);
% one Doppler
subplot(3,2,3)
plot(real(sigs(:, 2)));
title("1 Doppler (real)");
grid on;
xlim([1, N*fft_size]);
xticks(xindices);
subplot(3,2,4)
plot(imag(sigs(:, 2)));
title("1 Doppler (imag)");
grid on;
xlim([1, N*fft_size]);
xticks(xindices);
% one Delay
subplot(3,2,5)
plot(real(sigs(:, 3)));
title("1 delay (real)");
grid on;
xlim([1, N*fft_size]);
xticks(xindices);
subplot(3,2,6)
plot(imag(sigs(:, 3)));
title("1 delay (imag)");
grid on;
xlim([1, N*fft_size]);
xticks(xindices);

% low resolution
figure("Name", "Waveform (low resolution)")
xindices = 1:M:(N*M+1);
% base
subplot(3,2,1)
plot(real(sigs_base(:, 1)));
title("Base (real)");
grid on;
xlim([1, N*M]);
xticks(xindices);
subplot(3,2,2)
plot(imag(sigs_base(:, 1)));
title("Base (imag)");
grid on;
xlim([1, N*M]);
xticks(xindices);
% one Doppler
subplot(3,2,3)
plot(real(sigs_base(:, 2)));
title("1 Doppler (real)");
grid on;
xlim([1, N*M]);
xticks(xindices);
subplot(3,2,4)
plot(imag(sigs_base(:, 2)));
title("1 Doppler (imag)");
grid on;
xlim([1, N*M]);
xticks(xindices);
% one Delay
subplot(3,2,5)
plot(real(sigs_base(:, 3)));
title("1 delay (real)");
grid on;
xlim([1, N*M]);
xticks(xindices);
subplot(3,2,6)
plot(imag(sigs_base(:, 3)));
title("1 delay (imag)");
grid on;
xlim([1, N*M]);
xticks(xindices);


