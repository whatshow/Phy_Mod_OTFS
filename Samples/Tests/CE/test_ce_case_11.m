% Case 11:  All pilots allocation for center pilots
%           (Rectangular pulse & integer Doppler & reduced guard)
%
% Please note that there is no noise used in this example.
%
%% settings
test_ce_config;

%% all after channel
rgs = cell(length(Ms), 1);
rgs_rx = cell(length(Ms), 1);
for set_id = 1:length(Ms)
    % calculate the parameters
    N = Ns(set_id);
    M = Ms(set_id);
    pilots_num_delay = pilots_num_delays(set_id);
    pilots_num_doppl = pilots_num_doppls(set_id);
    %
    symbols_len = N*M-(pilots_num_delay+guard_delay_num_neg+guard_delay_num_pos)*(pilots_num_doppl+guard_doppl_num_neg+guard_doppl_num_pos);
    nbits = randi([0 1], symbols_len*M_bits, 1);
    xDD = qammod(nbits, M_mod,'InputType','bit','UnitAveragePower',true);
    % build rg
    rg = OTFSResGrid(M, N);
    rg.map(xDD, "pilots_num_delay", pilots_num_delay, "pilots_num_doppl", pilots_num_doppl, "pilots_pow", pil_pow, "guard_delay_num_neg", guard_delay_num_neg, "guard_delay_num_pos", guard_delay_num_pos, "guard_doppl_num_neg", guard_doppl_num_neg, "guard_doppl_num_pos", guard_doppl_num_pos);
    rg.setPulse2Recta();
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
    assert(sum(yDD_ful_diff)/N*M < 1e-10);
    % yDD diff
    H_DD = otfs.getChannel();
    yDD_diff = abs(yDD - H_DD*xDD);
    fprintf(" - Data RG yDD_diff: %e\n", sum(yDD_diff));
    assert(sum(yDD_diff)/symbols_len < 1e-10);
    % compare ce
    if pilots_num_delay == 1 && pilots_num_doppl == 1
        disp(" - the CE difference");
        fprintf("   - gain: %e\n", sum(abs(his_est - his)));
        fprintf("   - delay: %e\n", sum(lis-lis_est));
        fprintf("   - Doppler: %e\n", sum(kis-kis_est));
    end
end
%% plot
figure("name", "data & pilots & guards");
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