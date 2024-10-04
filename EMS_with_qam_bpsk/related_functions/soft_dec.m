function [HD_bin, HD1, yy] = soft_dec(p,N, y_cmp,norm_const,gray_bin_words_mod, I, Q )
HD1 = zeros(1, N);
HD_bin = zeros(N,p);
yy = HD_bin;

for j = 1 : N
    [HD_bin_mod, soft_dec, d, cndt, HD] = R_QAM_soft_dec(y_cmp(j), norm_const, gray_bin_words_mod, I, Q);
    HD_bin(j,:) = -(HD_bin_mod-1)/2;
    HD1(j) = bi2de(HD_bin(j,:),p);
    yy(j,:) = soft_dec;
end
