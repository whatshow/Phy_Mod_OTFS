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
N = 9;                          % time slot number (Doppler)
M = 11;                         % subcarrier number (Delay)
tap_pos = [5, 6; 9, 1; 1, 11];
delays = [1;0;1];
Doppler = [1;1;0];
descriptions = ["One Tap";"Doppler Overflow"; "Delay Overflow"];

% Generante graphs
for id = 1:length(delays)
    % get configuration
    tap_pos_x = tap_pos(id, 1);
    tap_pos_y = tap_pos(id, 2);
    li = delays(id);
    ki = Doppler(id);
    figure_title = descriptions(id);
    figure_title = figure_title + " (li = " + li + ", ki = " + ki + ")";
    
    % Gen information symbols (as a column vector)
    x_origin_DD = zeros(N, M);
    x_origin_DD(tap_pos_x, tap_pos_y) = sqrt(1/2)*(1 + 1j);
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
    % plot
    figure("Name", figure_title)
    subplot(1,2,1)
    bar3(abs(x_origin_DD));
    title("Before the channel")
    subplot(1,2,2)
    bar3(abs(Y_DD));
    title("After the channel")
end

