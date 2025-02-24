function u = code_decod_perm(x, h,add_mat, div_mat,per0)
N=length(x);
n=log2(N);
per0_inv(per0+1)=0:size(add_mat,1)-1;
m=nan(N,n+1);
m(:,end)=x;
l=n;
for t=0:N/2-1
    a=2*t-mod(t,2^(l-1))+1;
    b=2^(l-1)+2*t-mod(t,2^(l-1))+1;
    tmp1 = [m(a,l+1) m(b,l+1)];
    tmp1(2)=per0_inv(tmp1(2)+1);
    m(a,l)=add_mat(tmp1(1)+1, tmp1(2)+1);
    m(b,l)=tmp1(2);
end

for l=n-1:-1:1
    for t=0:N/2-1
        a=2*t-mod(t,2^(l-1))+1;
        b=2^(l-1)+2*t-mod(t,2^(l-1))+1;
         tmp1 = [m(a,l+1) m(b,l+1)];
         tmp1(2)=div_mat(tmp1(2)+1, h(b,l)+1);
         m(a,l)=add_mat(tmp1(1)+1, tmp1(2)+1);
         m(b,l)=tmp1(2);
    end
end
u=m(:,1);