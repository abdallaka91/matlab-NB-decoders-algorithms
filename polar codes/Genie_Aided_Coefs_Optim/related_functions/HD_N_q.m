function y = HD_N_q(LLR)
N = size(LLR,1);
y = zeros(N,1);
for i = 1 : N
    [~,n]=min(LLR(i,:));
    y(i,1) = n-1;
end