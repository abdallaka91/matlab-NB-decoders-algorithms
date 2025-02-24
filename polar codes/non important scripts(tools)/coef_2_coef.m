function  [coefs1, I1] = coef_2_coef(coefs)
N=size(coefs,1)*2;
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
h1(I1+1,:)=h;

coefs1=h;
