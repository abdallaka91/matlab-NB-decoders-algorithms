function [Pr0, Pr1]= F1(L1, L2, add_mat,mul_mat, h0, u0)
q = length(L1);
Pr0=zeros(q,1);
Pr1=zeros(q,1);


for u=0:q-1
    s=0;
    for up = 0:q-1
        a=add_mat(1+up, 1+u);
        b=mul_mat(1+h0, 1+up);
        s=s+L1(a+1)*L2(b+1);
    end
    Pr0(u+1)=s;
    a=add_mat(1+u0,1+u);
    b=mul_mat(1+h0, 1+u);
    Pr1(u+1)=L1(a+1)*L2(b+1);
end
Pr1=Pr1/sum(Pr1);
end
