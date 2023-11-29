%% OTFS Settings
N = 7;                          % time slot number
M = 12;                         % subcarrier number
N_syms_perfram = N*M;           % number of symbols per frame
N_bits_perfram = N*M*M_bits;    % number of bits per frame

%% Channel Settings
p = 18;                                         % channel paths number is set in every script
kmax = 1;                                       % maximal Doppler index
lmax = M - 1;                                   % maximal delay index (alll values are possible)

%% UE speed (km/h)
UEspeed = 500;