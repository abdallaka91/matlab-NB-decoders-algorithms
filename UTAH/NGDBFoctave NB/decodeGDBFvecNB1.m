function [d, failed_mx, Sgfo, Sbo, E, seqb, seqgf, l]  = decodeGDBFvecNB1(p, str_cn_vn,mul_mat, add_mat,y,hl, h, N, T, w, theta, nsigma)
M=size(h,1);
Nb = length(y);
d = sign(y);
dgf = HD_bin_gf(d,p);
Sgf = decod_prod(dgf,h,str_cn_vn, mul_mat, add_mat);
Sb = double(Sgf==0);
l=0;
failed_mx = M-sum(Sb);
Sgfo = Sgf;
Sbo=Sb;
% hold off
while l<T
    Sb1 = Sb;
    WSH = Sb1*hl;
    WSHb = zeros(1,Nb);
    for j = 1 : N
        WSHb((j-1)*p+1:j*p)=WSH(j);
    end

    E1  = d.*y + w*WSHb;
    E = E1+ nsigma*randn(1,Nb);
    
%     plot(E)
%     hold on
%     pause(0.005)
if l<100
[amp,flipdx]=mink(E,randi([10 50]));
elseif l<200
        flipdx = find(E<-0.7);
else
    flipdx = find(E<0.2);
end
    d(flipdx) = -d(flipdx);
    dgf = HD_bin_gf(d,p);
    Sgf = decod_prod(dgf,h,str_cn_vn, mul_mat, add_mat);
    Sb = double(Sgf==0);
    failed=M-sum(Sb);
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