function [Pr0, Pr1]= F1_permutation(L1, L2, add_mat,word1_p, u0)
q = length(L1);
Pr0=zeros(q,1);
Pr1=zeros(q,1);
word1_p_inv(word1_p+1)= 0:length(word1_p)-1;

%u (Pr0)-----------+-------------------------a   (L1)
%                  |
%                  |
%up (Pr1)----------=----<-----------Per->----b   (L2)

for a=0:q-1
    for b = 0:q-1
        bi=word1_p_inv(1+b);
        u=add_mat(1+a, 1+bi);
        c=L1(a+1)*L2(b+1);
        Pr0(u+1)= Pr0(u+1)+c;
    end
    up=add_mat(1+u0,1+a);
    b1=word1_p(1+up);
    Pr1(up+1)=L1(a+1)*L2(b1+1);
end
Pr1=Pr1/sum(Pr1);
end