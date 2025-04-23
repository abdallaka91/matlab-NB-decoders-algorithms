function [L1,gf, HD_L]=cond_prob_ccsk_1(y,sigm,etaq, k)
% L = -LLR_BPSK_GFq_2D(y, sigm, alph_bin);
q = size(etaq,1);
N = size(y,2);
L  = LLR_CCSK(y, etaq, q, N, sigm^2);
L1 = zeros(k, N);
gf = L1;
for i=1:N
    [L1(:,i), gf(:,i)] = mink(L(:,i), k);
end
gf = gf-1;
P1 = exp(-L);
P=bsxfun(@rdivide,P1,sum(P1,1));
[~,HD_L] = max(P,[],1);
HD_L = HD_L-1;