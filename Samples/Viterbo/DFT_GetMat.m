function dftmat = DFT_GetMat(N)
Wn = exp(-1j*2*pi/N);
dftmat = 1/sqrt(N)*Wn.^((0:N-1)'*(0:N-1));
end