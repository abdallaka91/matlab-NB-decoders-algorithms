function [d, failed_mx, Sgfo, Sbo, E, seqb, seqgf, l]  = decodeGDBFvecNB2(p, str_cn_vn,dc,mul_mat, add_mat,y,hl, h, N, T, w, theta, nsigma)
pw = 2.^(0:p-1);
theta0=theta;
stk=0;
M=size(h,1);
Nb = length(y);
d = sign(y);
dgf = HD_bin_gf(d,p,N,pw);
Sgf = decod_prod(dgf,h,str_cn_vn,dc, mul_mat, add_mat);
Sb = double(Sgf==0);
l=0;
failed_mx = M-sum(Sb);
Sgfo = Sgf;
Sbo=Sb;
% ofst_gf = sum(hl,1);
ofst_gf = mean(sum(hl,1))*ones(1,N);
ofst_b = 0.5*repmat(ofst_gf,p,1);
ofst_b=(ofst_b(:))';
% hold off
thk=1;
while l<T
    Sb1 = Sb;
    WSH = Sb1*hl;
    WSHb = zeros(1,Nb);
    for j = 1 : N
        WSHb((j-1)*p+1:j*p)=WSH(j);
    end

    E1  = d.*y + w*WSHb-ofst_b-1;
    E = E1+ nsigma*randn(1,Nb);
    
%     plot(E)
%     xlim([0 250])
%     hold on
%     pause(0.005)
    [a,i1] = mink(E,10);
    i1(a>theta0)=[];
    flipdx = i1;
    if theta0~=theta
        theta0=theta;
    end


    d(flipdx) = -d(flipdx);
    dgf = HD_bin_gf(d,p,N,pw);
    Sgf = decod_prod(dgf,h,str_cn_vn,dc, mul_mat, add_mat);
    Sb = double(Sgf==0);
    failed=M-sum(Sb);
    if failed<failed_mx
        stk=0;
        failed_mx = failed;
        Sgfo = Sgf;
        Sbo=Sb;
        seqb = (1-d)/2;
        seqgf = dgf;
        if failed_mx==0
            break
        end
    elseif failed_mx==failed
        stk=stk+1;
        if stk>100
            theta0= 0.3*theta;
            stk=0;
        end
    end

    l=l+1;
end
