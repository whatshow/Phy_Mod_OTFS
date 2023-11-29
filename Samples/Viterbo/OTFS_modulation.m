%% OTFS Modulation(VERIFIED)
% INPUT
% @x: x is N*M
%
% OUTPUT
% @s_mat: is M*N
% @s: each M elements are a column of total N colums 
function s = OTFS_modulation(N,M,x)
%% OTFS Modulation: 1. ISFFT, 2. Heisenberg transform
X = fft(ifft(x).').'/sqrt(M/N); %%%ISFFT
s_mat = ifft(X.')*sqrt(M); % Heisenberg transform
s = s_mat(:);

% X = ifft(x).'/sqrt(M/N); %%%ISFFT
% s_mat = X*sqrt(M); % Heisenberg transform
% s = s_mat(:);


%% Debug
% X2 = inv(dftmtx(N))*x.'/sqrt(M/N); %%%ISFFT
% s_mat2 = X*sqrt(M); % Heisenberg transform
% s2 = s_mat2(:);
% diff = sum(abs(s2 - s));
% fprintf("DEBUG: the difference between s2 and s is %f\n", diff);

end