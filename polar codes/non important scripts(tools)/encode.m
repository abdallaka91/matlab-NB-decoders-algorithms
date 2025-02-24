function x = encode(u, coefs,add_mat, mul_mat)
N=length(u);

i=-1;
u1=u;
while 1
    i=i+1;
    n=2^i;
    N1=N/n;
    if N1==1
        break
    end
    ns = reshape(1:N, N1, n);
    ns1= reshape(1:N/2, N/2/n, n);
    n1=size(ns,2);
    for j=1:n1
        cf1 = coefs(ns1(:,j),i+1);
        i2 = ns(:,j);
        u2 = u1(i2);
        u1(i2)=comb_func(u2,cf1, mul_mat, add_mat );
    end
end
x=u1;


function x = comb_func(u,cf, mul_mat, add_mat )
N = length(u);
x = zeros(size(u));
for i = 1 : N/2
    x(i)=add_mat(1+u(i), 1+u(i+N/2));
    x(i+N/2)=mul_mat(1+u(i+N/2), 1+cf(i));
end