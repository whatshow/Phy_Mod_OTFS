clear;
clc;

% QAM configuration
M_mod = 4;                                                                  % size of constellation
M_bits = log2(M_mod);
SNR = 10; % dB

% OTFS configuration
N = 9;                          % time slot number (Doppler)
M = 11;                         % subcarrier number (Delay)
tap_pos_x = 5;
tap_pos_y = 6;
li = 1;
ki = 1;

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
otfs.addChannelPath(1+1j, li, ki);
H_DD = otfs.getChannel();
% pass the channel
otfs.passChannel(0);
% demodulate
yDD = otfs.demodulate();
Y_DD = reshape(yDD, M, N).';

x = x_origin_DD(tap_pos_x, tap_pos_y);
y = Y_DD(tap_pos_x + ki, tap_pos_y + li);
gains = otfs.getChannelGains();
h = gains(1);
li_ki_est_rad = angle(y/h/x);
if li_ki_est_rad > 0
    li_ki_est_rad = li_ki_est_rad - 2*pi;
end
li_ki_est = li_ki_est_rad/(-2*pi/(N*M));

h_est_corr = y/x/exp(-1j*2*pi*li_ki_est/N/M);
h_est = y./x./exp(-1j*2*pi*ki/N*li/M);
h_est2 = y/x/exp(-1j*2*pi*(ki/N-li/M));
h_est3 = y/x/exp(-1j*2*pi*(li/M-ki/N));
h_est4 = abs(y/x)*exp(-1j*2*pi*ki*li/(M*N));
h_est5 = y/x/exp(1j*2*pi*ki*li/(M*N));
h_est6 = y/x/exp(1j*2*pi*ki*(tap_pos_y-1)/(M*N));

y_est = h*x*exp(-1j*2*pi*ki/N*li/M);
y_est2 = h_est2*x*exp(-1j*2*pi*li*ki/N/M);
y_est3 = h_est3*x*exp(-1j*2*pi*li*ki/N/M);
y_est4 = h_est4*x*exp(-1j*2*pi*li*ki/N/M);
y_est5 = h_est5*x*exp(-1j*2*pi*li*ki/N/M);
y_est6 = h_est6*x*exp(1j*2*pi*ki*(tap_pos_y-1)/(M*N));


% calculate the residual
fprintf("The channel phase equation is 2*pi*ki*(l-li)/M/N");
residual = abs(sum(y_est6 - y, "all"));
fprintf("The residual is %.16f\n", abs(residual));

% figure(1)
% subplot(2,2,1);
% bar3(real(x_origin_DD));
% title("X(real)");
% subplot(2,2,2);
% bar3(real(Y_DD));
% title("Y(real)");
% subplot(2,2,3);
% bar3(imag(x_origin_DD));
% title("X(imag)");
% subplot(2,2,4);
% bar3(imag(Y_DD));
% title("Y(imag)");

figure(2)
subplot(1,2,1);
bar3(abs(x_origin_DD));
title("X");
subplot(1,2,2);
bar3(abs(Y_DD));
title("Y");


clear Y_DD yDD H_DD;