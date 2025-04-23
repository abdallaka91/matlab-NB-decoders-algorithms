function [Pr0, Pr1, llr00, llr11, gf00, gf11]= F1_llr(L0, L1,gf0, gf1, add_mat,mul_mat, h0, u0, ofst, q)

mx0 = -1;
llr0=inf(q,1);
nm = length(gf0);


for i = 1 : nm
    u = gf0(i);
    for j=1:nm
        up = gf1(j);
        a=u;
        b=mul_mat(1+h0, 1+up);
        temp_gf = add_mat(1+a, 1+b);
        temp_llr = L0(i)+L1(j);
        if llr0(temp_gf+1)>temp_llr
            llr0(temp_gf+1)=temp_llr;
            if temp_llr>mx0
                mx0=temp_llr;
            end
        end
    end
end
llr0(llr0>1e30) = mx0+ofst;

L00=(ofst+L0(end))*ones(q,1);
L11=(ofst+L1(end))*ones(q,1);

L00(gf0+1)=L0;
L11(gf1+1)=L1;

llr1=inf(q,1);

for u=0:q-1
    a=add_mat(1+u0,1+u);
    b=mul_mat(1+h0, 1+u);
    llr1(u+1)=L00(a+1)+L11(b+1);
end

Pr0 = exp(-llr0);
Pr0 = Pr0/sum(Pr0);

Pr1 = exp(-llr1);
Pr1 = Pr1/sum(Pr1);


[llr00, gf00] = mink(llr0, nm);
gf00 = gf00-1;

[llr11, gf11] = mink(llr1, nm);
llr11 = llr11-llr11(1);
gf11 = gf11-1;