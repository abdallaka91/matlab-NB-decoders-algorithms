function V = ECN2(q, nm, nc, I, U, add_mat,comp_ECN) %        [m1,n1] = min(srt_llr);

I_gf = I(1,:);
I_llr = I(2,:);
U_gf = U(1,:);
U_llr = U(2,:);
V = U;
Mllr = zeros(nm, nm);
Mgf = Mllr;

for i = 1 : nm
    for j = 1 : nm
        Mllr(i,j) = I_llr(i)+U_llr(j);
        Mgf(i,j) = add_mat(I_gf(i)+1, U_gf(j)+1);
    end
end

nn=1;
V_llr = -ones(1,nm);
V_gf = -ones(1,nm);
srt_llr = Mllr(:,1);
srt_gf = Mgf(:,1);
crdn_vh = [1:nm; ones(1,nm);];
crdns = [1;1];
n1 = 1;
nnn=0;
while nn<=nm && nnn<=nc
    nnn=nnn+1;
    if crdns(2)<=nm
        srt_llr(n1) = Mllr(crdns(1),crdns(2));
        srt_gf(n1) = Mgf(crdns(1),crdns(2));
        [m1,n1] = min(srt_llr);
        gf1 = srt_gf(n1);
        crdns = crdn_vh(:,n1);
        crdns(2) = crdns(2)+1;
    else
        srt_llr(n1) = [];
        srt_gf(n1) = [];
        crdn_vh(:,n1)=[];
        [m1,n1] = min(srt_llr);
        gf1 = srt_gf(n1);
        crdns = crdn_vh(:,n1);
        crdns(2) = crdns(2)+1;
    end
    %     if sum(gf1==V_gf(1:nn))==0
    if ~any(gf1 == V_gf(1:nn))
        V_llr(nn) = m1;
        V_gf(nn) = srt_gf(n1);
        nn=nn+1;
    end
    crdn_vh(:,n1) =crdns;
end

if nn<=nm
    idx2 = find(V_gf<0);
    l1 = numel(idx2);
    if l1>0
        kk = max(V_llr);
        kk = kk + comp_ECN;
    end
    i=0;
    j=1;
    while i<= q-1 && l1>0
        if ~any(i == V_gf)
            V_gf(idx2(j)) = i;
            V_llr(idx2(j)) = kk;
            j = j+1;
            l1 = l1-1;
        end

        i=i+1;
    end
end
V(1,:) = V_gf;
V(2,:) = V_llr;
end