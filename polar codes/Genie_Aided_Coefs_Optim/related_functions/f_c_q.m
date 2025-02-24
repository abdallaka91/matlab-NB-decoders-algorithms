function Root = f_c_q(add_mat, mul_mat, Leafs, gf_coef)
N = size(Leafs,2);
m = N/2;
Root = zeros(1,N);
v1 = Leafs(1:m);
v2 = Leafs(m+1:N);
for i = 1 : m
    Root(i) = add_mat(v1(i)+1, v2(i)+1);
    Root(m+i:end) = mul_mat(v2(i)+1, gf_coef(i)+1);
end
% 

