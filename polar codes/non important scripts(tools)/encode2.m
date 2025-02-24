function [x,x1, I1] = encode2(u, coefs,add_mat, mul_mat)
N=length(u);
n=log2(N);
I1=bi2de(fliplr(de2bi((0:N-1)', n)));
N=length(u);
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
h1=nan(N,n);
h1(I1+1,:)=h;
u1(I1+1)=u;
m=nan(N,n+1);
m(:,1)=u;
for l=n:-1:1
    for t=0:N/2-1
        a=2*t-mod(t,2^(l-1))+1;
        b=2^(l-1)+2*t-mod(t,2^(l-1))+1;
        tmp1 = [m(a,n-l+1) m(b,n-l+1)];
        m(a,n-l+2)=add_mat(tmp1(1)+1, tmp1(2)+1);
        tmp2=[m(b,n-l+1) h(b,n-l+1)];
        m(b,n-l+2)=mul_mat(1+tmp2(1), 1+tmp2(2));
    end
end
x=m(:,end);

m1=nan(N,n+1);
m1(:,1)=u1;
for l=n:-1:1
    for t=0:N/2-1
        a=2*t-mod(t,2^(n-l))+1;
        b=2^(n-l)+2*t-mod(t,2^(n-l))+1;
        tmp1 = [m1(a,n-l+1) m1(b,n-l+1)];
        m1(a,n-l+2)=add_mat(tmp1(1)+1, tmp1(2)+1);
        tmp2=[m1(b,n-l+1) h1(b,n-l+1)];
        m1(b,n-l+2)=mul_mat(1+tmp2(1), 1+tmp2(2));
    end
end
x1=m1(:,end);