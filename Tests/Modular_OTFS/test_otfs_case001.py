import numpy as np
import sys
from OTFS import OTFS

M = 3; # subcarrier
N = 2; # timeslot
p = 2;
lmax = 1;
kmax = 1;
batch_size = 1;
constel = [-0.7071-0.7071j, -0.7071+0.7071j, 0.7071-0.7071j, 0.7071+0.7071j];

# no batch
# X_DD
xDD_idx = np.random.randint(4, size=(M*N));
xDD = np.take(constel,xDD_idx);
X_DD = np.reshape(xDD, (N, M));
# channel
otfs = OTFS();
for chy in range(2):
    if chy == 0:
        otfs.setChannel(p, lmax, kmax);
    else:
        otfs.setChannel(p, lmax, kmax, force_frac=True);
    
    # channel - ideal
    otfs.setPulse2Ideal();
    otfs.modulate(X_DD);
    otfs.passChannel(0);
    H_DD_ideal = otfs.getChannel();
    Y_DD_ideal = otfs.demodulate();
    # channel - recta
    otfs.setPulse2Recta();
    otfs.modulate(X_DD);
    otfs.passChannel(0);
    H_DD_recta = otfs.getChannel();
    Y_DD_recta = otfs.demodulate();
    # Y_DD
    yDD_ideal = np.reshape(Y_DD_ideal, M*N);
    yDD_recta = np.reshape(Y_DD_recta, M*N);
    
    # check residual
    yDD_ideal_diff = yDD_ideal - H_DD_ideal@xDD;
    assert(np.sum(abs(yDD_ideal_diff))/M/N < 20*np.finfo(float).eps);
    yDD_recta_diff = yDD_recta - H_DD_recta@xDD;
    assert(np.sum(abs(yDD_recta_diff))/M/N < 20*np.finfo(float).eps);  
# batch 1
batch_size = 1;
xDD_idx = np.random.randint(4, size=(batch_size, M*N));
xDD = np.take(constel,xDD_idx);
xDD_mat = np.expand_dims(xDD, axis=-1);
X_DD = np.reshape(xDD, (batch_size, N, M));
otfs = OTFS(batch_size=batch_size);
for chy in range(2):
    if chy == 0:
        otfs.setChannel(p, lmax, kmax);
    else:
        otfs.setChannel(p, lmax, kmax, force_frac=True);
    
    # channel - ideal
    otfs.setPulse2Ideal();
    otfs.modulate(X_DD);
    otfs.passChannel(0);
    H_DD_ideal = otfs.getChannel();
    Y_DD_ideal = otfs.demodulate();
    # channel - recta
    otfs.setPulse2Recta();
    otfs.modulate(X_DD);
    otfs.passChannel(0);
    H_DD_recta = otfs.getChannel();
    Y_DD_recta = otfs.demodulate();
    # Y_DD
    yDD_ideal = np.reshape(Y_DD_ideal, (batch_size, M*N, 1));
    yDD_recta = np.reshape(Y_DD_recta, (batch_size, M*N, 1));
    
    # check residual
    yDD_ideal_diff = yDD_ideal - H_DD_ideal@xDD_mat;
    assert(np.sum(abs(yDD_ideal_diff))/M/N/batch_size < 20*np.finfo(float).eps);
    yDD_recta_diff = yDD_recta - H_DD_recta@xDD_mat;
    assert(np.sum(abs(yDD_recta_diff))/M/N/batch_size < 20*np.finfo(float).eps);  
# batch n
batch_size = 8;
xDD_idx = np.random.randint(4, size=(batch_size, M*N));
xDD = np.take(constel,xDD_idx);
xDD_mat = np.expand_dims(xDD, axis=-1);
X_DD = np.reshape(xDD, (batch_size, N, M));
otfs = OTFS(batch_size=batch_size);
for chy in range(2):
    if chy == 0:
        otfs.setChannel(p, lmax, kmax);
    else:
        otfs.setChannel(p, lmax, kmax, force_frac=True);
    
    # channel - ideal
    otfs.setPulse2Ideal();
    otfs.modulate(X_DD);
    otfs.passChannel(0);
    H_DD_ideal = otfs.getChannel();
    Y_DD_ideal = otfs.demodulate();
    # channel - recta
    otfs.setPulse2Recta();
    otfs.modulate(X_DD);
    otfs.passChannel(0);
    H_DD_recta = otfs.getChannel();
    Y_DD_recta = otfs.demodulate();
    # Y_DD
    yDD_ideal = np.reshape(Y_DD_ideal, (batch_size, M*N, 1));
    yDD_recta = np.reshape(Y_DD_recta, (batch_size, M*N, 1));
    
    # check residual
    yDD_ideal_diff = yDD_ideal - H_DD_ideal@xDD_mat;
    assert(np.sum(abs(yDD_ideal_diff))/M/N/batch_size < 20*np.finfo(float).eps);
    yDD_recta_diff = yDD_recta - H_DD_recta@xDD_mat;
    assert(np.sum(abs(yDD_recta_diff))/M/N/batch_size < 20*np.finfo(float).eps);  
