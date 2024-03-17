% OTFS configuration
pilot_pow = 10^(50/10);

% test center - reduced guard - even size - even pilot number
N = 8;
M = 8;
pilots_num_delay = 2;
pilots_num_doppler = 2;
otfs = OTFS(M, N, "pilot_type", OTFS.PILOT_SINGLE_SISO, "pilot_loc_type", OTFS.PILOT_LOC_CENTER, "GUARD_TYPE", OTFS.GUARD_REDUCED);
otfs.insertPilotsAndGuards(pilots_num_delay, pilots_num_doppler, "guard_delay_num_neg", 1, "guard_delay_num_pos", 1, "pilots_pow", pilot_pow, "guard_doppler_num_neg", 1, "guard_doppler_num_pos", 1);
X_DD1 = otfs.X_DD;

% test center - reduced guard - even size - odd pilot number
N = 8;
M = 8;
pilots_num_delay = 1;
pilots_num_doppler = 1;
otfs = OTFS(M, N, "pilot_type", OTFS.PILOT_SINGLE_SISO, "pilot_loc_type", OTFS.PILOT_LOC_CENTER, "GUARD_TYPE", OTFS.GUARD_REDUCED);
otfs.insertPilotsAndGuards(pilots_num_delay, pilots_num_doppler, "guard_delay_num_neg", 1, "guard_delay_num_pos", 1, "pilots_pow", pilot_pow, "guard_doppler_num_neg", 1, "guard_doppler_num_pos", 1);
X_DD2 = otfs.X_DD;

% test center - reduced guard - odd size - even pilot number
N = 7;
M = 7;
pilots_num_delay = 2;
pilots_num_doppler = 2;
otfs = OTFS(M, N, "pilot_type", OTFS.PILOT_SINGLE_SISO, "pilot_loc_type", OTFS.PILOT_LOC_CENTER, "GUARD_TYPE", OTFS.GUARD_REDUCED);
otfs.insertPilotsAndGuards(pilots_num_delay, pilots_num_doppler, "guard_delay_num_neg", 1, "guard_delay_num_pos", 1, "pilots_pow", pilot_pow, "guard_doppler_num_neg", 1, "guard_doppler_num_pos", 1);
X_DD3 = otfs.X_DD;

% test center - reduced guard - odd size - odd pilot number
N = 7;
M = 7;
pilots_num_delay = 1;
pilots_num_doppler = 1;
otfs = OTFS(M, N, "pilot_type", OTFS.PILOT_SINGLE_SISO, "pilot_loc_type", OTFS.PILOT_LOC_CENTER, "GUARD_TYPE", OTFS.GUARD_REDUCED);
otfs.insertPilotsAndGuards(pilots_num_delay, pilots_num_doppler, "guard_delay_num_neg", 1, "guard_delay_num_pos", 1, "pilots_pow", pilot_pow, "guard_doppler_num_neg", 1, "guard_doppler_num_pos", 1);
X_DD4 = otfs.X_DD;

% plot
figure;
subplot(2,2,1);
bar3(abs(X_DD1));
title("reduced guard - even size - even pilot number")
subplot(2,2,2);
bar3(abs(X_DD2));
title("reduced guard - even size - odd pilot number")
subplot(2,2,3);
bar3(abs(X_DD3));
title("reduced guard - odd size - even pilot number")
subplot(2,2,4);
bar3(abs(X_DD4));
title("reduced guard - odd size - odd pilot number")