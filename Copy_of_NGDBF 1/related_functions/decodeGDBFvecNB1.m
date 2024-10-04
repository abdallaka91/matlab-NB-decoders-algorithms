function [d, failed_mx, Sgfo, Sbo, E, seqb, seqgf, l]  = decodeGDBFvecNB1(p,pw, h_arr,str_cn_vn,...
    dc,dv, mul_mat, add_mat,y,hl, h, N, T, w, theta, mnk,nsigma,dev_hamm1,dev_hamm2, sngl_VN, max_max_l)
rfrsh = 0;
stk=0;

M=size(h,1);
Nb = length(y);

dv1 = dv;
TT = -3.25*[1 1 1 ];

dv1(dv==2)=TT(1);
dv1(dv==3)=TT(2);
dv1(dv==4)=TT(3);

th1 = repmat(dv1,p,1);
th1 = (th1(:))';
th11 = th1;
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
    l = l + 1;
    WSH = bsxfun(@minus,w*Sb*hl,0.5*dv+1);
    for j = 1 : N
        WSHb((j-1)*p+1:j*p)=WSH(j);%-ofst_b(j);
    end

    E  = d.*y + WSHb+ nsigma*randn(1,Nb);

    flipdx = false(1,Nb);

    [a,i1] = sort(E);%mink(E,mnk);
    for i = 1 : length(i1)
        if E(i1(i))<th11(i1(i))
            flipdx(i1(i)) = true;
        else
            0;
        end
    end

    if rfrsh
        th11=th1;
        rfrsh = 0;
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
            th11= 0.3*th1;
            stk=0;
            rfrsh = 1;
        end
    end
end



if failed_mx>0 && sngl_VN
    fld_idx = Sbo==0;
    fld_vns=str_cn_vn(fld_idx);
    combinedArray = [fld_vns{:}]; % Combine all arrays into one
    uniqueIntegers = unique(combinedArray); % Find unique integers
    counts = histc(combinedArray, uniqueIntegers); % Count occurrences of each unique integer
    [~, i1] = sort(counts,'descend');
    vn1 = uniqueIntegers(i1);

    for o=1:min(3, length(vn1))
        dgf = seqgf;
        arr1 = y((vn1(o)-1)*p+1:(vn1(o)*p));
        % if counts(i1(o))==dv(vn1(o))
        cnddt = Lst_cndtt(arr1, dev_hamm1,pw);
        % elseif counts(i1(o))>1
        % cnddt = Lst_cndtt(arr1, dev_hamm2,pw);
        % end

        if counts(i1(o))==dv(vn1(o))
            for i = 1:numel(cnddt)
                l=l+1;
                if l > max_max_l
                    break;
                end
                dgf(vn1(o))=cnddt(i);
                Sgf = decod_prod1(dgf,h_arr,M, mul_mat, add_mat);
                Sb = double(Sgf==0);
                failed=M-sum(Sb);
                if failed<=failed_mx-dv(vn1(o))
                    failed_mx = failed;
                    Sgfo = Sgf;
                    Sbo=Sb;
                    seqgf = dgf;
                    break
                end
            end
            if failed_mx==0
                break;
            end
        end
        if l > max_max_l
            break;
        end
    end
end
end

