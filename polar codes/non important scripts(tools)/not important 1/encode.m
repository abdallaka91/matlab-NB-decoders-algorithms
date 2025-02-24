function x = encode(u, coefs,add_mat, mul_mat)
N=length(u);
% n=log2(N);
% t=0:N/2-1;
% h = ones(N,N);
% for l=n:-1:1
%     o=2*t-mod(t, 2^(l-1))+1;
%     f=2^(l-1)+2*t-mod(t, 2^(l-1))+1;
%     c = ones(N,1);
%     c(f)=coefs(:,l-1);
%     %     h(:,l+1) =  decod_prod(c,h,str_cn_vn, mul_mat, add_mat)
%     for i = 1 : N/2
%         h(f(i),l-1) = mul_mat(c(f(i))+1, h(f(i),l)+1);
%     end
% end
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