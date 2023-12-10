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
sigs_X_DD = zeros(length(descriptions), N, M);
sigs_X_DT = zeros(length(descriptions), N*M);
for id = 1:length(descriptions)
    % get configuration
    tap_pos_x = tap_pos(id, 1);
    tap_pos_y = tap_pos(id, 2);
    
    X_DD = zeros(N, M);
    X_DD(tap_pos_x, tap_pos_y) = sqrt(1/2)*(1 + 1j);
    % store the data
    sigs_X_DD(id, :, :) = X_DD;
    
    % init OTFS
    otfs = OTFS(M, N);
    % modulate
    otfs.modulate(X_DD);
    % get the symbol in the time domain
    sigs(:, id) = otfs.getS("fft_size", fft_size);
    sigs_base(:, id) = otfs.getS("fft_size", M);
    X_DT = otfs.getXDT();
    X_DT = X_DT(:);
    sigs_X_DT(id, :) = X_DT;
end


%% plot
% in the delay Doppler domain
figure("Name", "Data in DD domain")
for id = 1:length(descriptions)
    subplot(1,length(descriptions),id);
    bar3(abs(squeeze(sigs_X_DD(id, :, :))));
    title(descriptions(id));
    ylabel("Doppler");
    xlabel("delay")
end
% time domain (high resolution)
figure("Name", "Waveform (high resolution)")
xindices = 1:fft_size:(N*fft_size+1);
for id = 1:length(descriptions)
    subplot(length(descriptions),2,2*id-1)
    plot(real(sigs(:, id)));
    title(descriptions(id) + " (real)");
    grid on;
    xlim([1, N*fft_size]);
    xticks(xindices);
    subplot(length(descriptions),2,2*id)
    plot(imag(sigs(:, id)));
    title(descriptions(id) + " (imag)");
    grid on;
    xlim([1, N*fft_size]);
    xticks(xindices);
end

% time domain (low resolution)
figure("Name", "Waveform (low resolution)")
xindices = 1:M:(N*M+1);
% base
for id = 1:length(descriptions)
    subplot(length(descriptions),2,2*id-1)
    plot(real(sigs_base(:, id)));
    title(descriptions(id) + " (real)");
    grid on;
    xlim([1, N*M]);
    xticks(xindices);
    subplot(length(descriptions),2,2*id)
    plot(imag(sigs_base(:, id)));
    title(descriptions(id) + " (imag)");
    grid on;
    xlim([1, N*M]);
    xticks(xindices);
end

% delay time domain (HD)
figure("Name", "Waveform in delay time domain (high resolution)")
xindices = 1:M:(N*M+1);
for id = 1:length(descriptions)
    subplot(length(descriptions),2,2*id-1)
    plot(real(sigs_X_DT(id, :)));
    title(descriptions(id) + " (real)");
    grid on;
    xlim([1, N*M]);
    xticks(xindices);
    subplot(length(descriptions),2,2*id)
    plot(imag(sigs_X_DT(id, :)));
    title(descriptions(id) + " (imag)");
    grid on;
    xlim([1, N*M]);
    xticks(xindices);
end


