%% OTFS Demodulation (VERIFIED)
% INPUT
% @r: each M elements are a column of total N colums
%
% OUTPUT
% @y: is a N*M matrix
function y = OTFS_demodulation(N,M,r)
%% OTFS demodulation: 1. Wiegner transform, 2. SFFT
%% Original Viterbo method
r_mat = reshape(r,M,N);
Y_MN = fft(r_mat)/sqrt(M); % Wigner transform
Y_NM = Y_MN.';
y = ifft(fft(Y_NM).').'/sqrt(N/M); % SFFT

%% Clearer Viterbo method
% r_mat = reshape(r,M,N);
% Y = fft(r_mat)/sqrt(M); % Wigner transform
% y = ifft(fft(Y, [], 2)).'/sqrt(N/M); % SFFT


%% Debug (always on M*N)
% r_mat_debug = reshape(r,M,N);
% Y_debug = fft(r_mat_debug)/sqrt(M); % Wigner transform
% y_debug1 = fft(ifft(Y_debug), [], 2).'/sqrt(N/M); % SFFT
% 
% y_debug2 = fft(r_mat_debug, [], 2).'/sqrt(N/M)/sqrt(M);
% y_debug3 = (r_mat_debug*dftmtx(N)).'/sqrt(N/M)/sqrt(M);
% 
% diff_y = abs(y_debug1 - y);
% diff_y_debug_1_2 = abs(y_debug1 - y_debug2);
% diff_y_debug_2_3 = abs(y_debug3 - y_debug2);
% 
% fprintf("DEBUG: done!\n");
end