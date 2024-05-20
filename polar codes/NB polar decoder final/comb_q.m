function y = comb_q(v1, v2, add_mat)
N1 = length(v1);
y = zeros(1, 2*N1);
for i = 1 : N1
    temp = add_mat(v1(i)+1, v2(i)+1);
    y([i,i+N1]) =[temp, v2(i)];
end