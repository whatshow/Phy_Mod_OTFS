%% OTFS Delay-Dopple Channel Generation (Effective)
% This channel is used to offer the calculate effective channel for the
% vectorized input and output
% 
% INPUT
% @M: 
% @N:
% @taps: the number of delay doppler channel taps
% @delay_taps:      those taps' coordinate in the delay axis
% @Doppler_taps:    those taps' coordinate in the doppler axis
% @chan_coef:       those taps' coefficient
% @sigma_2:         the channel variance in the time-frequency domain
%
% OUTPUT
% @H:           the time-frequency channel
% @Heff_rect:   the delay-doppler channel for rectangular pulses
function [H, Heff_rect] = OTFS_ddchannel_gen(M, N, taps,delay_taps,doppler_taps,chan_coef, varargin)
    % Inputs Name-Value Pair 
    inPar = inputParser;
    addParameter(inPar, 'real_channel', false, @islogical);
    inPar.KeepUnmatched = true;         % Allow unmatched cases
    inPar.CaseSensitive = false;        % Allow capital or small characters
    parse(inPar, varargin{:});
    real_channel = inPar.Results.real_channel;

    % Output
    H = zeros(M*N, M*N);
    Heff_rect = zeros(M*N, M*N);

    % DFT matrix and IDFT matrix
    
    dftmat = dftmtx(N)/sqrt(N);
    idftmat = conj(dftmat);
       
    % Permutation Matrix
    piMat = eye(M*N);
    for itap = 1:taps
        li_val = delay_taps(itap);
        ki_val = doppler_taps(itap);
        
        hi = chan_coef(itap);
        %piMati = piMat^li_val;  % Viterbo's way to get the permutation matrix
        piMati = circshift(piMat, li_val);  % My way to get the permutation matrix
        
        ki = ki_val;
        ki_idx = circshift(-li_val:M*N-1-li_val, -li_val);
        deltaMat_diag = exp(1j*2*pi*ki/(M*N)*ki_idx);
        deltaMati = diag(deltaMat_diag);
        %deltaMati = double(subs(deltaMat, ki, ki_val));
        
        deltaMati2 = diag(exp(1j*2*pi*ki/(M*N)*(-li_val:M*N-1-li_val)));
        
        % Calculate the channel in the time-frequency domain
        %H = H + hi*piMati*deltaMati;
        H = H + hi*deltaMati2*piMati;
        %H(abs(H) < 1e-13) = 0;
        
        % Calculate the channel in the delay doppler domain
        % Calculate Pi(delay) and Qi(doppler)
        Pi = kron(dftmat, eye(M))*piMati*kron(idftmat, eye(M));
        %Pi = kron(dftmat, eye(M))*piMati;
        %Qi = kron(dftmat, eye(M))*deltaMati*kron(idftmat, eye(M));
        Qi2 = kron(dftmat, eye(M))*deltaMati2*kron(idftmat, eye(M));
        %Qi = deltaMati*kron(idftmat, eye(M));
        
        % Calculate Ti
        %Ti = Pi*Qi;
        Ti = Qi2*Pi;
        
        % update Heff_rect
        Heff_rect = Heff_rect + hi*Ti;
    end
end