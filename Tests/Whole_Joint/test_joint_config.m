% OTFS
fc = 4;                 % GHz
freq_sp = 15;           % kHz
M = 32;                 % time slot
N = 16;                 % subcarrier
% speed
vs = [30, 120, 500];    % km/h
fd = vs/3.6/physconst('LightSpeed')*fc*1e6; % Doppler shift (kHz)
kmaxs = fd/(freq_sp/N);
lmaxs = 4*ones(length(kmaxs), 1);
% SNR range
SNR_ds = 10:2:18;       % dB
SNR_ps = 25:5:50;       % dB
% pilots
pl_len = 1;
pk_len = 1;
gdn_len = ceil(max(lmaxs));
gdp_len = ceil(max(lmaxs));
% QAM
M_mod = 4;
M_bits = log2(M_mod);
constel = qammod(0: M_mod-1, M_mod, 'UnitAveragePower',true);
% OTFS frame
N_syms_perfram = N*M-(pl_len + gdn_len + gdp_len)*N;
N_bits_perfram = N_syms_perfram*M_bits;
% Monte Carlo
N_Frams = 1e6*ones(length(SNR_ds), 1);
% file
path_folder = "./_data/Samples/Tests/Whole_Joint/";
% create the folder if not exist
if ~exist(path_folder, 'dir')
    mkdir(path_folder);
end

