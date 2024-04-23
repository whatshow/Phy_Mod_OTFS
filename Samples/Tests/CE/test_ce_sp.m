% Superimposed pilots
%
%% settings
test_ce_config;

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
    symbols_len = N*M;
    nbits = randi([0 1], symbols_len*M_bits, 1);
    xDD = qammod(nbits, M_mod,'InputType','bit','UnitAveragePower',true);
    % build rg
    rg = OTFSResGrid(M, N);
    rg.setPulse2Ideal();
    rg.setPilot2SuperImposed();
    rg.setPilot2Center(pilots_num_delay, pilots_num_doppl);
    rg.setGuard(guard_delay_num_neg, guard_delay_num_pos, guard_doppl_num_neg, guard_doppl_num_pos);
    rg.map(xDD, "pilots_pow", pil_pow);
    rgs{set_id} = rg;
    % through the channel
    otfs = OTFS();
    otfs.modulate(rg);
    otfs.setChannel(p, lmax, kmax);
    otfs.passChannel(0);
    rg_rx = otfs.demodulate();
    rgs_rx{set_id} = rg_rx;
end
%% plot
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