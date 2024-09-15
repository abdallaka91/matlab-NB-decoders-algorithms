function [failed_mx, seqgf,seqgf_cmplx, l]  = decodeGDSFvecNB_QAM(h_arr,mul_mat, add_mat,...
    y,hl, h, T,max_max_l, str_vn_cn, str_cn_vn,dc, dv, w, theta,nsigma,mnk,  cnstl1, code_seq)

M=size(h,1);
N = size(h,2);
[HD0_gf, HD0_cplx] = HD_QAM(y, cnstl1);
Synd0_gf = decod_prod1(HD0_gf',h_arr,M, mul_mat, add_mat);
Synd0_bin = double(Synd0_gf==0);
failed_mx = M-sum(Synd0_bin);

l=0;

if failed_mx==0
    seqgf = HD0_gf;
else
    Synd_bin = Synd0_bin;
    HD_cplx = HD0_cplx;
    HD_gf = HD0_gf;
    while l<T
        l = l + 1;
        nse =  nsigma*(randn(N,1)+1i*randn(N,1));
        dst = abs(nse);
        y1 = y + nse;
        [HD_nse_gf, HD_nse_cplx]  = HD_QAM(y1, cnstl1);
        WSH = Synd_bin*hl;
        E = WSH - w*dst.';
        %         plot(E,'.')
        i2 = find(E<theta);
        [~,i3] = mink(E(i2),mnk);
        i1=i2(i3);
        HD_cplx(i1) = HD_nse_cplx(i1);
        HD_gf(i1) = HD_nse_gf(i1);

        Synd_gf = decod_prod1(HD_gf',h_arr,M, mul_mat, add_mat);
        Synd_bin = double(Synd_gf==0);
        failed = M-sum(Synd_bin);
        if failed < failed_mx
            Synd_bin0 = Synd_bin;
            seqgf_cmplx = HD_cplx;
            seqgf = HD_gf;
            failed_mx = failed;
            if failed==0
                break
            end
        end

    end
    if failed_mx>10000
        fld_idx = find(Synd_bin0==0);
        fld_vns=str_cn_vn(fld_idx);
        combinedArray = [fld_vns{:}]; % Combine all arrays into one
        uniqueIntegers = unique(combinedArray); % Find unique integers
        counts = histc(combinedArray, uniqueIntegers); % Count occurrences of each unique integer
        [~, i1] = sort(counts,'descend');
        vn1 = uniqueIntegers(i1);

        for o=1:min(3, length(vn1))
            dgf = seqgf;
            dgf_cmplx = seqgf_cmplx;
            [~,i2] = sort(abs(y(o)-cnstl1));
            cnddt_cmplx = cnstl1(i2);
            cnddt_gf = i2-1;
            if counts(i1(o))==dv(vn1(o))
                for i = 1:numel(cnddt_gf)
                    l=l+1;
                    if l > max_max_l
                        break;
                    end
                    dgf_cmplx(vn1(o))=cnddt_cmplx(i);
                    dgf(vn1(o))=cnddt_gf(i);

                    Synd_gf = decod_prod1(dgf',h_arr,M, mul_mat, add_mat);
                    Synd_bin = double(Synd_gf==0);
                    failed=M-sum(Synd_bin);
                    if failed<=failed_mx-dv(vn1(o))+1
                        failed_mx = failed;
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




