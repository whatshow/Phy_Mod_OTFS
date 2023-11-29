%% OTFS QAM Modulation
%
%
function x = OTFS_QamMod(bits,M_mod, isUnitaryPower)
    if isUnitaryPower == true
        x = qammod(bits,M_mod,'gray', 'UnitAveragePower', true);
    else
        x = qammod(bits,M_mod,'gray');
    end
end