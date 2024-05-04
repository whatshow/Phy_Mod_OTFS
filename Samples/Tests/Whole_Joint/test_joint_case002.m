clear;
clc;
test_joint_config;
%% rewrite configuration
% channel
p = 9;
lmax = lmaxs(1);
kmax = kmaxs(end);
% data powers
SNR_d = SNR_ds(end);
pow_sig = 1;                                % signal power
pow_pils = 10.^((SNR_ps - SNR_d)/10);       % pilot power
No = 10.^(-SNR_d/10);                       % noise power
pow_thr = 3*sqrt(No);                       % threshold
N_Frams = 1e4*ones(length(pow_pils), 1);
path_file = path_folder + "test_joint_case001" + ".mat";

%% simulations
SER = zeros(length(SNR_ps), 1);
SER2 = zeros(length(SNR_ps), 1);
for idx = 1:length(SNR_ps)
    % retrieve info
    SNR_p = SNR_ps(idx);
    pow_pil = pow_pils(idx);
    N_Fram = N_Frams(idx);
    % print
    fprintf("Pilot SNR is %d\n", SNR_p);
    %
    tmp_SER = zeros(N_Fram, 1);
    tmp_SER2 = zeros(N_Fram, 1);
    parfor i_Fram = 1:N_Fram
        % build RG
        data_temp = randi([0 1], N_bits_perfram, 1);
        xDD = qammod(data_temp, M_mod, 'InputType', 'bit', 'UnitAveragePower', true);
        rg = OTFSResGrid(M, N);
        rg.setPilot2Center(pl_len, pk_len);
        rg.setGuard(gdn_len, gdp_len, "guard_doppl_full", true);
        rg.setPulse2Recta();
        rg.map(xDD, "pilots_pow", pow_pil);
        rg2 = OTFSResGrid(M, N);
        rg2.setPilot2Center(pl_len, pk_len);
        rg2.setGuard(gdn_len, gdp_len, "guard_doppl_full", true);
        rg2.setPulse2Ideal();
        rg2.map(xDD, "pilots_pow", pow_pil);
        % through the channel
        otfs = OTFS();
        otfs.modulate(rg);
        otfs.setChannel(p, lmax, kmax);
        otfs.passChannel(No);
        [his, lis, kis] = otfs.getCSI("sort_by_delay_doppler", true);   % perfect CSI CSI
        rg_rx = otfs.demodulate();
        otfs.modulate(rg2);
        otfs.passChannel(No);
        rg2_rx = otfs.demodulate();
        % demapping (CE)
        [yDD, his_est, lis_est, kis_est] = rg_rx.demap("isData", false, "threshold", pow_thr);
        [yDD2, his_est2, lis_est2, kis_est2] = rg2_rx.demap("isData", false, "threshold", pow_thr);
        % compare the channel
        HDD = otfs.getChannel(his, lis, kis);
        HDD_norm2 = sum(abs(HDD).^2, "all");
        if isempty(his_est)
            tmp_H_NMSE(i_Fram) = 1;
        else
            HDD_est = otfs.getChannel(his_est, lis_est, kis_est);
            tmp_H_NMSE(i_Fram) = ;
        end
        if isempty(his_est2)
            tmp_H_NMSE2(i_Fram) = 10*log10(HDD_norm2);
        else
            HDD_est2 = otfs.getChannel(his_est2, lis_est2, kis_est2);
            tmp_H_NMSE2(i_Fram) = 10*log10(sum(abs(HDD-HDD_est2).^2, "all")/HDD_norm2);
        end
    end
    H_NMSE(idx) = mean(tmp_H_NMSE);
    H_NMSE2(idx) = mean(tmp_H_NMSE2);
end
%% save
save(path_file);

%% plot
semilogy(SNR_ps, H_NMSE, "-s", "Color", "#D95319", "LineWidth", 4);
hold on;
semilogy(SNR_ps, H_NMSE2, "--ob", "LineWidth", 2);
hold off;
grid on;
xlabel("Pilot SNR(dB)");
ylabel("NMSE(dB)");
ylim([min(H_NMSE2), max(H_NMSE)]);
xlim([min(SNR_ps), max(SNR_ps)]);
legend('rect', 'ideal');
title(string(M)+"x"+string(N)+", "+string(p)+" paths (Full Guard)")