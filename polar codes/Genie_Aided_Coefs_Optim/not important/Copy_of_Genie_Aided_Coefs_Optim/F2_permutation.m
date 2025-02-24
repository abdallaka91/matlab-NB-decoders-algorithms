function [Pr0, Pr1]= F2_permutation(L1, L2, add_mat,per0,per1, u0)
q = length(L1);
Pr0=zeros(q,1);
Pr1=zeros(q,1);
per0_inv(per0+1)= 0:length(per0)-1;
per1_inv(per1+1)= 0:length(per1)-1;

%u (Pr0)-----------+----------------Per0->---------a   (L1)
%                  |
%                  |
%up (Pr1)----------=----------------Per1->---------b   (L2)

for a=0:q-1
    for b = 0:q-1
        bi=per1_inv(1+b);
        ai=per0_inv(1+a);
        u=add_mat(1+ai, 1+bi);
        c=L1(a+1)*L2(b+1);
        Pr0(u+1)= Pr0(u+1)+c;
    end
    ai=per0_inv(1+a);
    up=add_mat(1+u0,1+ai);
    b1=per1(1+up);
    Pr1(up+1)=L1(a+1)*L2(b1+1);
end
Pr1=Pr1/sum(Pr1);
end