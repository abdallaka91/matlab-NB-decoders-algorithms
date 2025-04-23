function [u,x,m1, y]=gen_noisy_seq_ccsk(N,q,sigm, add_mat, etaqm,u)
if nargin==5
    u = randi([0 q-1], 1, N);
end
m1=encode_layer_n_ccsk(u,add_mat);
x=m1(:,end)';
etax = zeros(q, N);
for i = 1 : N
    etax(:,i)=etaqm(:,x(i)+1);
end

nse = sigm*randn(size(etax));
y = etax + nse;
