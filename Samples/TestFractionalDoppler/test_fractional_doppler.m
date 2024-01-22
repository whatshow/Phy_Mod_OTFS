clear;
clc;

% QAM configuration
M_mod = 4;                                                                  % size of constellation
M_bits = log2(M_mod);
sympool = qammod(0: M_mod-1, M_mod, 'UnitAveragePower',true);               % Generate the symbol pool
sympool_real = unique(real(sympool));
SNR = 10; % dB
No = 0;

% OTFS configuration
N = 9;                          % time slot number
M = 11;                         % subcarrier number
p = 1;
lmax = 2;
DopplerValue = [2, 3.6];
description = ["Integer Doppler", "Fractional Doppler"];

%% Gen information symbols (as a column vector)
sigs = zeros(length(description), N, M);
sigs_X_DD = zeros(length(description), N, M);
sigs_dopplers = [0, 0];
for id = 1:length(description)
    x_origin_DD = sqrt(1/100)*(ones(N, M) + 1j*ones(N, M));
    x_origin_DD(:, 6-2:6+2) = 0;
    x_origin_DD(5, 6) = sqrt(1/2)*(1 + 1j);
    x_origin = x_origin_DD.';
    x_origin = x_origin(:);
    
    % init OTFS
    otfs = OTFS(M, N);
    % modulate
    otfs.modulate(x_origin_DD);
    % set the channel
    H_DD = otfs.setChannel("p", p, "kmax", DopplerValue(id), "lmax", lmax);
    dopplers = otfs.getChannelDopplers();
    sigs_dopplers(id) = dopplers;
    % pass the channel
    otfs.passChannel(No);
    % demodulate
    yDD = otfs.demodulate();
    Y_DD = reshape(yDD, M, N).';
    
    % calculate the residual
    residual = sum(yDD - H_DD*x_origin, "all");
    fprintf("The residual is %.16f\n", abs(residual));
    
    % store
    sigs(id, :, :) = Y_DD;
    sigs_X_DD(id, :, :) = x_origin_DD;
end

%% plot
for id = 1:length(description)
    % plot
    figure("Name", description(id) + " kmax=" + num2str(sigs_dopplers(id)))
    subplot(1,2,1)
    bar3(abs(squeeze(sigs_X_DD(id, :, :))));
    title("Before the channel")
    subplot(1,2,2)
    bar3(abs(squeeze(sigs(id, :, :))));
    title("After the channel")
end


