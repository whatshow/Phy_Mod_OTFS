clear;
clc;

M = 3; % subcarrier
N = 2; % timeslot
p = 2;
lmax = 1;
kmax = 1;

% X_DD
constel = [-0.7071-0.7071j, -0.7071+0.7071j, 0.7071-0.7071j, 0.7071+0.7071j].';
xDD_idx = randi(4, M*N, 1);
xDD = constel(xDD_idx);
X_DD = reshape(xDD, M, N).';
% channel
otfs = OTFS();
for chy = 1:2
    if chy == 1
        otfs.setChannel(p, lmax, kmax);
    else
        otfs.setChannel(p, lmax, kmax, "force_frac", true);
    end
    
    % channel - ideal
    otfs.setPulse2Ideal();
    otfs.modulate(X_DD);
    otfs.passChannel(0);
    H_DD_ideal = otfs.getChannel();
    Y_DD_ideal = otfs.demodulate();
    % channel - recta
    otfs.setPulse2Recta();
    otfs.modulate(X_DD);
    otfs.passChannel(0);
    H_DD_recta = otfs.getChannel();
    Y_DD_recta = otfs.demodulate();
    % Y_DD
    yDD_ideal = reshape(Y_DD_ideal.', M*N, 1);
    yDD_recta = reshape(Y_DD_recta.', M*N, 1);
    
    % check residual
    yDD_ideal_diff = yDD_ideal - H_DD_ideal*xDD;
    assert(sum(abs(yDD_ideal_diff))/M/N < 20*eps);
    yDD_recta_diff = yDD_recta - H_DD_recta*xDD;
    assert(sum(abs(yDD_recta_diff))/M/N < 20*eps);
end