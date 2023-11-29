clear;
clc;

% QAM configuration
M_mod = 4;                                                                  % size of constellation
M_bits = log2(M_mod);
sympool = qammod(0: M_mod-1, M_mod, 'UnitAveragePower',true);               % Generate the symbol pool
sympool_real = unique(real(sympool));

% OTFS configuration
N = 7;                          % time slot number
M = 12;                         % subcarrier number
N_syms_perfram = N*M;           % number of symbols per frame
N_bits_perfram = N*M*M_bits;    % number of bits per frame


% Gen information symbols (as a column vector)
data_info_bit = randi([0,1],N_bits_perfram,1);
data_temp = bi2de(reshape(data_info_bit,N_syms_perfram,M_bits));
x_origin = qammod(data_temp,M_mod,'gray', 'UnitAveragePower', true);

% init OTFS
otfs = OTFS(M, N);