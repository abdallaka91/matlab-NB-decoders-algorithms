function [d, failed_mx, Sgfo, Sbo, E, seqb, seqgf, l]  = decodeGDBFvecNB4(p,pw, h_arr,str_cn_vn,dc,dv, mul_mat, add_mat,y,hl, h, N, T, w, theta, mnk,nsigma)
theta0=theta;
stk=0;
M=size(h,1);
Nb = length(y);
d = sign(y);
seqb = 0.5*(1-d);
dgf = HD_bin_gf(d,p,N,pw);
seqgf= dgf;
Sgf = decod_prod1(dgf,h_arr,M, mul_mat, add_mat);
Sb = double(Sgf==0);
l=0;
failed_mx = M-sum(Sb);
Sgfo = Sgf;
Sbo=Sb;

WSHb=zeros(1,Nb);
while l<T

    WSH = w*Sb*hl-0.5*mean(dv)-1;
    for j = 1 : N
        WSHb((j-1)*p+1:j*p)=WSH(j);%-ofst_b(j);
    end

    E  = d.*y + WSHb+ nsigma*randn(1,Nb);

    [a,i1] = mink(E,mnk);
    i1(a>theta0)=[];
    flipdx = i1;
    if theta0~=theta
        theta0=theta;
    end


    d(flipdx) = -d(flipdx);
    dgf = HD_bin_gf(d,p,N,pw);
    Sgf = decod_prod1(dgf,h_arr,M, mul_mat, add_mat);
    Sb = double(Sgf==0);
    failed=M-sum(Sb);
    if failed<failed_mx || failed==0
        do =d;
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

%--------

if failed_mx>0
    fld_idx = find(Sbo==0);
    fld_vns=str_cn_vn(fld_idx);
    combinedArray = [fld_vns{:}]; % Combine all arrays into one
    uniqueIntegers = unique(combinedArray); % Find unique integers
    counts = histc(combinedArray, uniqueIntegers); % Count occurrences of each unique integer
    [~, i1] = sort(counts,'descend');
    vn1 = uniqueIntegers(i1);
    % dv1 = dv(vn1);
    % ls1 = nan(1,dv1);
    % c1=0;
    % for i = 1 : length(fld_idx)
    %     if sum(vn1==fld_vns{i})>0
    %         c1=c1+1;
    %         ls1(c1)=fld_idx(i);
    %
    %     end
    % end
    % if c1~=dv1
    %     ls1(isnan(ls1))=[];
    % end

    for o=1:3
        dgf = seqgf;
        if counts(i1(o))==dv(vn1(o))
            for cnddt=1:2^p
                l=l+1;
                dgf(vn1(o))=cnddt-1;
                Sgf = decod_prod1(dgf,h_arr,M, mul_mat, add_mat);
                Sb = double(Sgf==0);
                failed=M-sum(Sb);
                if failed<failed_mx
                    failed_mx = failed;
                    Sgfo = Sgf;
                    Sbo=Sb;
                    seqgf = dgf;
                    if failed_mx==0
                        break
                    end
                end
            end
            if failed_mx==0
                break;
            end
        end
    end
end
end

