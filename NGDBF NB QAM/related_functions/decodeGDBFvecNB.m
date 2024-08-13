function [d, failed_mx, Sgfo, Sbo, E, seqb, seqgf, l]  = decodeGDBFvecNB(p,pw, h_arr,str_cn_vn,...
    dc,dv, mul_mat, add_mat,y,hl, h, N, T, w, theta, mnk,nsigma,dev_hamm1,dev_hamm2, sngl_VN,...
    max_max_l, gray_labels, cod_seqq)
theta0=theta;
stk=0;
M=size(h,1);
Nb = length(y);
d = sign(y);
seqb = 0.5*(1-d);
dgfff = HD_bin_gf(d,p,N,pw);
dgf = gray_labels(dgfff+1);
seqgf= dgf;
Sgf = decod_prod1(dgf,h_arr,M, mul_mat, add_mat);
Sb = double(Sgf==0);
l=0;
failed_mx = M-sum(Sb);
Sgfo = Sgf;
Sbo=Sb;

if failed_mx==0
    E=[];
else
    WSHb=zeros(1,Nb);
    while l<T
        l = l + 1;
        WSH = bsxfun(@minus,w*Sb*hl,0.5*dv+1);
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
        Sgf = decod_prod1(gray_labels(dgf+1),h_arr,M, mul_mat, add_mat);
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
                    Sgf = decod_prod1(gray_labels(dgf+1),h_arr,M, mul_mat, add_mat);
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
end

