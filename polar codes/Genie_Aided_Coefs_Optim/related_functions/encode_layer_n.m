function m1=encode_layer_n_ccsk(u,add_mat, mul_mat,l1)
N=length(u);
n=log2(N);

m1=nan(N,n+1);
m1(:,l1)=u;
for l=l1:n
    for t=0:N/2-1
        a=2*t-mod(t,2^(l-1))+1;
        b=2^(l-1)+2*t-mod(t,2^(l-1))+1;
        tmp1 = [m1(a,l-1+1) m1(b,l-1+1)];
        m1(a,l-1+2)=add_mat(tmp1(1)+1, tmp1(2)+1);
        m1(b,l-1+2)=m1(b,l-1+1);
    end
end
