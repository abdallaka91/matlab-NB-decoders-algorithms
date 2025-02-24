function [P, HD_L]=cond_prob(y,sigm,alph_bin)
L = -LLR_BPSK_GFq_2D(y, sigm, alph_bin);
P1 = exp(-L);
P=bsxfun(@rdivide,P1,sum(P1,1));
[~,HD_L] = max(P,[],1);
HD_L = HD_L-1;