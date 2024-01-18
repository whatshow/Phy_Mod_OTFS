import numpy as np
from numpy.random import default_rng
import matplotlib.pyplot as plt
import sys 
sys.path.append("..") 
from OTFS import OTFS;

# QAM configuration
sympool = np.asarray([-0.707106781186548+0.707106781186548j, -0.707106781186548-0.707106781186548j,0.707106781186548+0.707106781186548j,0.707106781186548-0.707106781186548j]);
SNR = 10; # dB
#No = 1/10^(SNR/10); # linear
No = 0;

# OTFS configuration
N = 9;                          # time slot number
M = 11;                         # subcarrier number
delays = [1, 1];
Doppler = [1, 1.43];
description = ["Full Guards(Integer Doppler)", "Full Guards(Fractional Doppler)"];


## Gen information symbols (as a column vector)
sigs = np.zeros((len(description), N, M), dtype=np.complex128);
sigs_X_DD = np.zeros((len(description), N, M), dtype=np.complex128);


for des_id in range(len(description)):
    hi = 1;
    li = delays[des_id];
    ki = Doppler[des_id];
    
    x_origin_DD = np.sqrt(1/100)*(np.ones((N, M)) + 1j*np.ones((N, M)));
    x_origin_DD[:, 5-2:5+2+1] = 0;
    x_origin_DD[4, 5] = np.sqrt(1/2)*(1 + 1j);
    x_origin = np.reshape(x_origin_DD, N*M);
    
    # init OTFS
    otfs = OTFS(M, N);
    # modulate
    otfs.modulate(x_origin_DD);
    # set the channel
    otfs.addChannelPath(hi, li, ki);
    H_DD = otfs.getChannel();
    # pass the channel
    otfs.passChannel(No);
    # demodulate
    yDD = otfs.demodulate();
    Y_DD = np.reshape(yDD, (N, M));
    
    # calculate the residual
    residual = abs(yDD - np.matmul(H_DD, x_origin));
    print("The residual is %.16f\n"%np.sum(residual, axis=None));
    
    # store
    sigs[des_id, :, :] = Y_DD;
    sigs_X_DD[des_id, :, :] = x_origin_DD;
    
# plot
for des_id in range(len(description)):
    x = [];
    y = [];
    z = [];
    dz = [];
    dz2 = [];
    for N_id in range(N):
        for M_id in range(M):
            x.append(N_id - 0.5);
            y.append(M_id - 0.5);
            z.append(0);
            dz.append(abs(sigs_X_DD[des_id, N_id, M_id]));
            dz2.append(abs(sigs[des_id, N_id, M_id]));
    dx = 0.5*np.ones_like(x);
    dy = 0.5*np.ones_like(y);
    
    fig = plt.figure(des_id, figsize=(8, 3), constrained_layout=True); # constrained_layout avoid overlapping titles
    fig.suptitle(description[des_id]);
    ax1 = fig.add_subplot(121, projection='3d');
    ax2 = fig.add_subplot(122, projection='3d');
    ax1.bar3d(x, y, z, dx, dy, dz, shade=True);
    ax1.set_title('Before the channel')
    ax2.bar3d(x, y, z, dx, dy, dz2, shade=True);
    ax2.set_title('After the channel')
    plt.show();
    