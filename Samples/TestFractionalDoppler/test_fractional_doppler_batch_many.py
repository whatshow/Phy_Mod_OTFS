import numpy as np
from numpy.random import default_rng
import matplotlib.pyplot as plt
import sys 
sys.path.append("..") 
from OTFS import OTFS;

# QAM configuration
sympool = np.asarray([-0.707106781186548+0.707106781186548j, -0.707106781186548-0.707106781186548j,0.707106781186548+0.707106781186548j,0.707106781186548-0.707106781186548j]);
SNR = 10; # dB
No = 0;


# OTFS configuration
N = 9;                          # time slot number
M = 11;                         # subcarrier number
p = 1;
lmax = 2;
DopplerTag = ["kmax", "kmax_frac"];
DopplerValue = [2, 3.6];
description = ["Integer Doppler", "Fractional Doppler"];

# batch
batch_size = 2;

# Gen information symbols (as a column vector)
sigs =  np.zeros((len(description), batch_size, N, M), dtype=np.complex128);
sigs_X_DD = np.zeros((len(description), batch_size, N, M), dtype=np.complex128);
sigs_dopplers = [0, 0];
for des_id in range(len(description)):
    x_origin_DD = np.sqrt(1/100)*(np.ones((N, M)) + 1j*np.ones((N, M)));
    x_origin_DD[:, 5-2:5+2+1] = 0;
    x_origin_DD[4, 5] = np.sqrt(1/2)*(1 + 1j);
    x_origin_DD = np.tile(x_origin_DD, (batch_size, 1, 1));
    x_origin = np.reshape(x_origin_DD, (batch_size, N*M));
    
    # init OTFS
    otfs = OTFS(M, N, batch_size=batch_size);
    # modulate
    otfs.modulate(x_origin_DD);
    # set the channel
    otfs.setChannel(p=p, kmax=DopplerValue[des_id], lmax=lmax);
    H_DD = otfs.getChannel();
    dopplers = otfs.getChannelDopplers();
    sigs_dopplers[des_id] = dopplers;
    # pass the channel
    otfs.passChannel(No);
    # demodulate
    yDD = otfs.demodulate();
    Y_DD = np.reshape(yDD, (batch_size, N, M));
    
    # calculate the residual
    residual = abs(yDD - np.squeeze(np.matmul(H_DD, np.expand_dims(x_origin, axis=-1)), axis=-1));
    print("The residual is %.16f\n"%np.sum(residual, axis=None));
    
    # store
    sigs[des_id, ...] = Y_DD;
    sigs_X_DD[des_id, ...] = x_origin_DD;

# plot
for des_id in range(len(description)):
    fig = plt.figure(des_id, figsize=(8, (1+8*batch_size)), constrained_layout=True); # constrained_layout avoid overlapping titles
    figure_size = batch_size*100+22;
    fig.suptitle(description[des_id]);
    for batch_id in range(batch_size):
        cur_figure_id = (batch_id + 1)*100;
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
                dz.append(abs(sigs_X_DD[des_id, batch_id, N_id, M_id]));
                dz2.append(abs(sigs[des_id, batch_id, N_id, M_id]));
        dx = 0.5*np.ones_like(x);
        dy = 0.5*np.ones_like(y);
        ax1 = fig.add_subplot((cur_figure_id + 21), projection='3d');
        ax2 = fig.add_subplot((cur_figure_id + 22), projection='3d');
        ax1.bar3d(x, y, z, dx, dy, dz, shade=True);
        ax1.set_title('Before the channel')
        ax2.bar3d(x, y, z, dx, dy, dz2, shade=True);
        ax2.set_title('After the channel, kmax=' + str(sigs_dopplers[des_id][batch_id][0]))
    plt.show();