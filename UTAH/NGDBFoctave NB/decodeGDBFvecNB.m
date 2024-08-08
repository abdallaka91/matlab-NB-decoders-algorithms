function [d, failed_mx, Sgfo, Sbo, E, seqb, seqgf, l]  = decodeGDBFvecNB(p, str_cn_vn,mul_mat, add_mat,y,hl, h, N, T, w, theta, nsigma)
Nb = length(y);
d = sign(y);
dgf = HD_bin_gf(d,p);
Sgf = decod_prod(dgf,h,str_cn_vn, mul_mat, add_mat);
Sb = double(Sgf~=0);
l=0;
failed_mx = sum(Sb);
Sgfo = Sgf;
Sbo=Sb;
while l<T
    Sb1 = 1-2*Sb;
    WSH = Sb1*hl;
    WSHb = zeros(1,Nb);
    for j = 1 : N
        WSHb((j-1)*p+1:j*p)=WSH(j);
    end

    E = d.*y + w*WSHb + nsigma*randn(1,Nb);
    flipdx = find(E<theta);
    d(flipdx) = -d(flipdx);
    dgf = HD_bin_gf(d,p);
    Sgf = decod_prod(dgf,h,str_cn_vn, mul_mat, add_mat);
    Sb = double(Sgf~=0);
    failed=sum(Sb);
    if failed<failed_mx
        failed_mx = failed;
        Sgfo = Sgf;
        Sbo=Sb;
        seqb = (1-d)/2;
        seqgf = dgf;
        if failed_mx==0
            break
        end
    end
    l=l+1;
end
