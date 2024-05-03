% Case 22:  All pilots allocation for center pilots
%           (Ideal pulse & fractional Doppler & full Guard)
%
%% settings
test_ce_config;
kmax = 1.8;

%% pass the channel
rgs = cell(length(Ms), 1);
rgs_rx = cell(length(Ms), 1);
for set_id = 1:length(Ms)
    % calculate the parameters
    N = Ns(set_id);
    M = Ms(set_id);
    pilots_num_delay = pilots_num_delays(set_id);
    pilots_num_doppl = pilots_num_doppls(set_id);
    % 
    symbols_len = N*M-(pilots_num_delay+guard_delay_num_neg+guard_delay_num_pos)*N;
    nbits = randi([0 1], symbols_len*M_bits, 1);
    xDD = qammod(nbits, M_mod,'InputType','bit','UnitAveragePower',true);
    % build rg
    rg = OTFSResGrid(M, N);
    rg.setPulse2Ideal();
    
    rg.setPilot2Center(pilots_num_delay, pilots_num_doppl);
    rg.setGuard(guard_delay_num_neg, guard_delay_num_pos, 0, 0,"guard_doppl_full", true);
    rg.map(xDD, "pilots_pow", pil_pow);
    rgs{set_id} = rg;
    % through the channel
    otfs = OTFS();
    otfs.modulate(rg);
    otfs.setChannel(p, lmax, kmax);
    otfs.passChannel(0);
    rg_rx = otfs.demodulate();
    rgs_rx{set_id} = rg_rx;
    [yDD, his_est, lis_est, kis_est] = rg_rx.demap("threshold", 1e-10);
    [his, lis, kis] = otfs.getCSI("sort_by_delay_doppler", true);
    %% data check
    fprintf("No.%0d: OTFS %dx%d, Pilot %dx%d\n", set_id, N, M, pilots_num_doppl, pilots_num_delay);
    % HDD (full)
    H_DD_ful_est = otfs.getChannel(his, lis, kis, "data_only", false);
    H_DD_ful = otfs.getChannel("data_only", false);
    H_DD_ful_diff = abs(H_DD_ful - H_DD_ful_est);
    fprintf(" - HDD (full) diff: %e\n", sum(H_DD_ful_diff, "all"));
    assert(sum(H_DD_ful_diff, "all") < eps);
    % yDD (full) diff
    xDD_ful = rg.getContent("isVector", true);
    yDD_ful = rg_rx.getContent("isVector", true);
    yDD_ful_diff = abs(yDD_ful - H_DD_ful*xDD_ful);
    fprintf(" - Full RG yDD_diff: %e\n", sum(yDD_ful_diff));
    assert(sum(yDD_ful_diff)/N*M < eps);
    % yDD diff
    H_DD = otfs.getChannel();
    yDD_diff = abs(yDD - H_DD*xDD);
    fprintf(" - Data RG yDD_diff: %e\n", sum(yDD_diff));
    assert(sum(yDD_diff)/symbols_len < 1e-13);
    % compare channel
    if pilots_num_delay == 1 && pilots_num_doppl == 1
        est_H_DD_ful = otfs.getChannel(his_est, lis_est, kis_est, "data_only", false);
        est_H_DD = otfs.getChannel(his_est, lis_est, kis_est);
        est_H_DD_ful_diff = abs(est_H_DD_ful - H_DD_ful);
        est_H_DD_diff = abs(est_H_DD - H_DD);
        [est_H_DD_diff_r, est_H_DD_diff_c] = size(est_H_DD_diff);
        disp(" - Estimated Channel");
        fprintf("   - HDD (full) diff %e\n", sum(est_H_DD_ful_diff, "all")/M/N/M/N);
        assert(sum(est_H_DD_ful_diff, "all")/M/N/M/N < eps);
        fprintf("   - HDD diff %e\n", sum(est_H_DD_diff, "all")/est_H_DD_diff_r/est_H_DD_diff_c/est_H_DD_diff_r/est_H_DD_diff_c);
        assert(sum(est_H_DD_diff, "all")/est_H_DD_diff_r/est_H_DD_diff_c < eps);
    end
end
% plot
% plot - before the channel
fig_title = "RG";
figure("name", fig_title);
for set_id = 1:length(Ms)
    % calculate the parameters
    N = Ns(set_id);
    M = Ms(set_id);
    pilots_num_delay = pilots_num_delays(set_id);
    pilots_num_doppl = pilots_num_doppls(set_id);
    % plot
    subplot(2,2,set_id);
    bar3(abs(rgs{set_id}.getContent()));
    ylabel("Doppler");
    xlabel("delay");
    title("OTFS:" + string(N) + "x" + string(M) + ", Pilot:" + string(pilots_num_doppl) + "x" + string(pilots_num_delay));
end
sgtitle(fig_title);
% plot - after the channel
fig_title = "RG (after channel)";
figure("name", fig_title);
for set_id = 1:length(Ms)
    % calculate the parameters
    N = Ns(set_id);
    M = Ms(set_id);
    pilots_num_delay = pilots_num_delays(set_id);
    pilots_num_doppl = pilots_num_doppls(set_id);
    % plot
    subplot(2,2,set_id);
    bar3(abs(rgs_rx{set_id}.getContent()));
    title("OTFS:" + string(N) + "x" + string(M) + ", Pilot:" + string(pilots_num_doppl) + "x" + string(pilots_num_delay));
end
sgtitle(fig_title);