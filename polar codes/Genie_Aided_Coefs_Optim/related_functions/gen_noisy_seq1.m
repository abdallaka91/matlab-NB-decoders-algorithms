function [u,x,m1, y]=gen_noisy_seq1(N,q,sigm,perms_n, add_mat, alph_bin_mod,l1, u)
if nargin==7
    u = randi([0 q-1], 1, N);
end
 m1=encode_layer_n_permut1(u,add_mat,perms_n, l1);
 x=m1(:,end)';
etax = alph_bin_mod(x+1, :);
nse = sigm*randn(size(etax));
y = etax + nse;