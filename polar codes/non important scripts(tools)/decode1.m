function u = decode1(x, coefs,add_mat, div_mat)
N=length(x);
n=log2(N);
I1=bi2de(fliplr(de2bi((0:N-1)', n)));
h=nan(N,n);
for l=1:n
    b=nan(1,N/2);
    for t=0:N/2-1
        a=2*t-mod(t,2^(l-1));
        b(t+1)=2^(l-1)+2*t-mod(t,2^(l-1));
    end
    b1=reshape(b, 2^(l-1),N/2/(2^(l-1)))';
    b2=b1(:);
    h(b2+1,n-l+1)=coefs(:,n-l+1);
end

m=nan(N,n+1);
m(:,end)=x;
for l=n:-1:1
    for t=0:N/2-1
        a=2*t-mod(t,2^(n-l))+1;
        b=2^(n-l)+2*t-mod(t,2^(n-l))+1;
         tmp1 = [m(a,l+1) m(b,l+1)];
         tmp1(2)=div_mat(tmp1(2)+1, h(b,l)+1);
         m(a,l)=add_mat(tmp1(1)+1, tmp1(2)+1);
         m(b,l)=tmp1(2);
    end
end
u=m(:,1);