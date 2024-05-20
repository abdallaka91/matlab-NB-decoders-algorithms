function Root = f_c_q(add_mat, Leafs)
N = size(Leafs,2);
m = N/2;
Root = zeros(1,N);
v1 = Leafs(1:m);
v2 = Leafs(m+1:N);
for i = 1 : m
    Root(i) = add_mat(v1(i)+1, v2(i)+1);
end
Root(m+1:end) = v2;

