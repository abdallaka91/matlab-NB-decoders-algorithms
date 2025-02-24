function m1=encode_layer_n_permut1(u,add_mat,word1_pc, l1)
N=length(u);
n=log2(N);

m1=nan(N,n+1);
m1(:,l1)=u;
for l=l1:n
    word1_pl=word1_pc{l};
    bb=0:N/(2^l)-1;
    bbb=repmat(bb, (N/2)/length(bb),1);
    bbb=bbb(:);
    for t=0:N/2-1
%         b=mod(t,N/(2^l))+1;
        word1_p=word1_pl(bbb(t+1)+1,:);
        a=2*t-mod(t,2^(l-1))+1;
        b=2^(l-1)+2*t-mod(t,2^(l-1))+1;
        tmp1 = [m1(a,l-1+1) m1(b,l-1+1)];
        m1(a,l-1+2)=add_mat(tmp1(1)+1, tmp1(2)+1);
%         tmp2=[m1(b,l-1+1) h1(b,l-1+1)];
%         m1(b,l-1+2)=mul_mat(1+tmp2(1), 1+tmp2(2));
        m1(b,l-1+2)=word1_p(1+tmp1(2));
    end
end
