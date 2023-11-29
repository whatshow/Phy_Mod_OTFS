%% OTFS QAM Demodulation
%
%
function bits = OTFS_QamDemod(syms,M_mod, isUnitaryPower)
    if isUnitaryPower == true
        bits = qamdemod(syms,M_mod,'gray', 'UnitAveragePower', true);
    else
        bits = qamdemod(syms,M_mod,'gray');
    end
    
end