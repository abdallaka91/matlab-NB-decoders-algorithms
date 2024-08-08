function V = ECN4(q, nm, I, U, add_mat) %        [m1,n1] = min(srt_llr);

I_gf = I(1,:)+1;
I_llr = I(2,:);
U_gf = U(1,:)+1;
U_llr = U(2,:);
V = U;
Mllr = zeros(nm*nm,1);
Mgf = Mllr;
k=0;
for i = 1 : nm
    for j = 1 : nm
        k = k+1;
        Mllr(k) = I_llr(i)+U_llr(j);
        Mgf(k) = add_mat(I_gf(i), U_gf(j))+1;
    end
end

V_gf = false(1,q);
V_llr = zeros(1,q);

[Mllr1,b] = sort(Mllr);
Mgf1 = Mgf(b);

k=1;
for i = 1 : nm^2
    if  (V_gf(Mgf1(i)))==false
        V_gf(Mgf1(i)) = true;
        V_llr(Mgf1(i)) = Mllr1(i);
        k=k+1;
    end
    if k> nm
        break;     
    end
end
V_gf1=find(V_gf);
V_llr1 = V_llr(V_gf1);
[V_llr1,b1] = sort(V_llr1);
V_gf = V_gf1(b1)-1;
V(1,:) = V_gf;
V(2,:) = V_llr1(1:nm);
end
