clear;
clc;
c_len=4;

c = (1:c_len)*(1+1j).';
c = [1 0 0 0];
C = zeros(c_len, c_len);
for col_id = 1:c_len
    C(:,col_id)=circshift(c, col_id-1);
end

F = dftmtx(c_len);
Fn = conj(F)/c_len;
U = F/sqrt(c_len);
Un = Fn*sqrt(c_len);

C_diag_1 = U*C*Un;
C_diag_1_2 = Fn*C*F;
C_diag_2 = diag(fft(c));
residual_C_diag_1 = C_diag_1 - C_diag_1_2;
residual_C_diag_12 = C_diag_1 - C_diag_2;

residual_c = fft(c.')/sqrt(c_len) - F*(c.');