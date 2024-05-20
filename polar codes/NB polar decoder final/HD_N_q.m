function y = HD_N_q(prb)
N = size(prb,1);
y = zeros(N,1);
for i = 1 : N
    [~,n]=max(prb(i,:));
    y(i,1) = n-1;
end