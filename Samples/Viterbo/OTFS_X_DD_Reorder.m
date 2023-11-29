% reorder x_DD to get the binary sequence
function x_DD_ordered = OTFS_X_DD_Reorder(M, N, x_DD)
    x_DD_MN = reshape(x_DD,M,N);
    x_DD_NM = x_DD_MN.';
    x_DD_ordered = x_DD_NM(:);
end