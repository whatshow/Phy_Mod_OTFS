clear;
clc;
% load settings
test_mp_config;

%% simulation
conste = qammod(0:M_mod-1,M_mod, 'gray');
parfor i_fram = 1:10000
    % random input bits generation
    data_info_bit = randi([0,1],N_bits_perfram,1);
    data_temp = bi2de(reshape(data_info_bit,N_syms_perfram,M_bits));
    x = qammod(data_temp,M_mod,'gray');
    x = reshape(x,N,M);
    % OTFS modulation
    s = OTFS_modulation(N,M,x);
    % OTFS channel generation
    taps = 4;
    chan_coef = [-0.7323+0.5284i,0.0746-0.1501i,-0.2060+0.5961i,0.0064+0.1292i];
    delay_taps = [0, 1, 1, 1];
    Doppler_taps = [1, 0, -1, -1];
    % OTFS channel output
    L = max(delay_taps);
    s = [s(N*M-L+1:N*M);s];%add one cp
    s_chan = 0;
    for p_id = 1:taps
        s_chan = s_chan+chan_coef(p_id)*circshift([s.*exp(1j*2*pi/M ...
            *(-L:-L+length(s)-1)*Doppler_taps(p_id)/N).';zeros(delay_taps(end),1)],delay_taps(p_id));
    end
    r = s_chan + noise;
    r = r(L+1:L+(N*M));%discard cp
    
    % OTFS demodulation
    y = OTFS_demodulation(N,M,r);
    
    %% detection
    x_est = OTFS_mp_detector(N,M,M_mod,taps,delay_taps,Doppler_taps,chan_coef,sigma_2,y);
    x_est2 = OTFS_MP(N,M, taps,chan_coef, delay_taps,Doppler_taps,sigma_2,y, conste);

    x_est_diff = abs(x_est2 - x_est);
    assert(sum(x_est_diff, "all") == 0);
end
disp("MP code is the same.");