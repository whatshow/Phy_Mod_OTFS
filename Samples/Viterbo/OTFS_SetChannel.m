%% Set the Channel
% <INPUTS>
% Pp: the number of paths
% @lmax: the maximal delay
% @kmax: the maximal doppler
function [delay_taps,Doppler_taps,chan_coef] = OTFS_SetChannel(P, lmax, kmax)

% the minimal delay is 1
lmin= 1;

% the minima Doppler is the negative of maximal Doppler
kmin = - kmax;

% generate all possible paths when we only select P paths among them and
% give them channel gains (others are 0)
taps_max = (kmax - kmin + 1)*(lmax - lmin + 1);
delay_taps_all = kron(lmin:lmax, ones(1, kmax - kmin + 1)); % create [lmin, lmin, lmin, lmin+1, lmin+1, lmin+1 ...]
Doppler_taps_all = repmat(kmin:kmax, 1, lmax - lmin + 1);

% We only suppose there are P paths

% We select P paths from all possible paths; that is, we do the
% randperm(taps_max) and we choose the first P items
taps_idx_chaotic = randperm(taps_max);
taps_selected_idx = taps_idx_chaotic(1:P);


pow_prof = (1/P) * (ones(1, P));
chan_coef = sqrt(pow_prof).*(sqrt(1/2) * (randn(1, P)+1i*randn(1, P)));

% Then we evalue the varable to the output
delay_taps = delay_taps_all(taps_selected_idx);
% set minimal delay = 0
delay_taps(find(delay_taps == min(delay_taps), 1)) = 0;
Doppler_taps = Doppler_taps_all(taps_selected_idx);

end