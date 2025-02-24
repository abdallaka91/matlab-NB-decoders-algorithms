function x1 = encode_coef_perm(u, h1,add_mat, mul_mat, perm)
N=length(u);
n=log2(N);

m1=nan(N,n+1);
m1(:,1)=u;
for l=1:n-1
    for t=0:N/2-1
        a=2*t-mod(t,2^(l-1))+1;
        b=2^(l-1)+2*t-mod(t,2^(l-1))+1;
        tmp1 = [m1(a,l-1+1) m1(b,l-1+1)];
        m1(a,l-1+2)=add_mat(tmp1(1)+1, tmp1(2)+1);
        tmp2=[m1(b,l-1+1) h1(b,l-1+1)];
        m1(b,l-1+2)=mul_mat(1+tmp2(1), 1+tmp2(2));
    end
end
l=n;
for t=0:N/2-1
    a=2*t-mod(t,2^(l-1))+1;
    b=2^(l-1)+2*t-mod(t,2^(l-1))+1;
    tmp1 = [m1(a,l-1+1) m1(b,l-1+1)];
    m1(a,l-1+2)=add_mat(tmp1(1)+1, tmp1(2)+1);
    m1(b,l-1+2)=perm(1+m1(b,l-1+1));
end

x1=m1(:,end);