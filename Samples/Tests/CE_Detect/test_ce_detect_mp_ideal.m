clear;
clc;

%% File
path_folder = "./_data/Samples/Tests/CE_Detect/";
path_file = path_folder + "test_ce_detect_mp_ideal.mat";
% create the folder if not exist
if ~exist(path_folder, 'dir')
    mkdir(path_folder);
end
% delete the file if exist
if exist(path_file, 'file')
    delete(path_file)
end


%% settings
SNR_p = 30;                         % Pilots SNR (dB)
SNR_ds = 4:2:20;                     % Data SNR (dB)
Nos = 10.^(-SNR_ds/10);               % noise power
pil_pows = 10.^((SNR_p - SNR_ds)/10);
%pil_pow = Nos.*10.^(SNR_p/10);
pil_thrs = 3*sqrt(Nos);             % pilot threshold
% OTFS configuration
N = 12;                             % time slot number
M = 90;                             % subcarrier number
% channel
p = 4;
lmax = 3;
kmax = 6;
% pilot settings
pilots_num_delay = 1;
pilots_num_doppl = 1;
guard_delay_num_neg = lmax;
guard_delay_num_pos = lmax;
% QAM configuration
M_mod = 4;
M_bits = log2(M_mod);
sympool = qammod(0: M_mod-1, M_mod, 'UnitAveragePower',true);
% OTFS frame
N_syms_perfram = N*M-(pilots_num_delay+guard_delay_num_neg+guard_delay_num_pos)*N;
N_bits_perfram = N_syms_perfram*M_bits;
% Monte Carlo
N_Frams = 3e4*ones(length(SNR_ds), 1);
N_Frams(end-4:end-3)=1e5;
N_Frams(end-2)=5e5;
N_Frams(end-1:end)=1e6;

%% simulations
SERs_PerCSI = zeros(length(SNR_ds), 1);
SERs_CE = zeros(length(SNR_ds), 1);
for SNR_d_id = 1:length(SNR_ds)
    SNR_d = SNR_ds(SNR_d_id);
    pil_thr = pil_thrs(SNR_d_id);
    pil_pow = pil_pows(SNR_d_id);
    No = Nos(SNR_d_id);
    fprintf("SNR(data) is %d\n", SNR_d);
    N_Fram = N_Frams(SNR_d_id);
    tmp_SERs_PerCSI = zeros(N_Fram, 1);
    tmpSERs_CE = zeros(N_Fram, 1);
    parfor i_Fram = N_Fram
        % build RG
        data_temp = randi([0 1], N_bits_perfram, 1);
        xDD = qammod(data_temp, M_mod, 'InputType', 'bit', 'UnitAveragePower', true);
        rg = OTFSResGrid(M, N);
        rg.setPulse2Ideal();
        rg.map(xDD, "pilots_num_delay", pilots_num_delay, "pilots_num_doppl", pilots_num_doppl, "pilots_pow", pil_pow, "guard_delay_num_neg", guard_delay_num_neg, "guard_delay_num_pos", guard_delay_num_pos, "guard_doppl_full", true);
        % through the channel
        otfs = OTFS();
        otfs.modulate(rg);
        otfs.setChannel(p, lmax, kmax);
        otfs.passChannel(No);
        [his, lis, kis] = otfs.getCSI("sort_by_delay_doppler", true);   % perfect CSI CSI 
        rg_rx = otfs.demodulate();
        % demapping (CE)
        [~, his_est, lis_est, kis_est] = rg_rx.demap("isData", false, "threshold", pil_thr);
        % detect
        od = OTFSDetector(sympool);
        od.useMPBase();
        xDD_est_percsi  = od.detect(rg_rx, his, lis, kis, No);
        xDD_est_ce      = od.detect(rg_rx, his_est, lis_est, kis_est, No);
        tmp_SERs_PerCSI(i_Fram) = sum(abs(xDD_est_percsi - xDD) > eps)/N_syms_perfram;
        if empty(his_est)
            tmpSERs_CE(i_Fram) = 1;
        else
            tmpSERs_CE(i_Fram) = sum(abs(xDD_est_ce - xDD) > eps)/N_syms_perfram;
        end
    end
    % update the average SER
    SERs_PerCSI(SNR_d_id) = mean(tmp_SERs_PerCSI);
    SERs_CE(SNR_d_id) = mean(tmpSERs_CE);
    % save
    save(path_file);
end
% save
save(path_file);

%% plot
semilogy(SNR_ds, SERs_CE, "-s", "Color", "#D95319", "LineWidth", 4);
hold on;
semilogy(SNR_ds, SERs_PerCSI, "--ob", "LineWidth", 2);
hold off;
grid on;
xlabel("SNR(dB)");
ylabel("SER");
ylim([min(SERs_alva), 1]);
xlim([min(SNR_range), max(SNR_range)]);
legend('PerCSI-MP', 'CE-MP');