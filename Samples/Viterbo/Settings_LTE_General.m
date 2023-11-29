%% QAM Settings
M_mod = 4;                                                                  % size of constellation
M_bits = log2(M_mod);
sympool = qammod(0: M_mod-1, M_mod, 'UnitAveragePower',true);               % Generate the symbol pool
sympool_real = [real(sympool(1)), imag(sympool(1))];

%% Simulation Settings
SNR_dB = 8;                                % SNR Range
SNR = 10.^(SNR_dB/10);
noise_var_sqrt = sqrt(1./SNR);
eng_sqrt = 1;                                   % average energy per data symbol
sigma_2 = abs(eng_sqrt*noise_var_sqrt).^2;      % noiser power

N_fram = 3e4*ones(length(SNR_dB), 1).';         % frame number to eliminate the burst error effect (more frames for 14dB, 16dB)

%% Detection Settings
iter_times = 10;    % BPIC, AMP or EP iteration times

%% EP beta
ep_beta = 0.9;
