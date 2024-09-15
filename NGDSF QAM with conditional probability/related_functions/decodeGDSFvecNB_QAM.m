function [failed_mx, seqgf, l]  = decodeGDSFvecNB_QAM(h_arr,mul_mat, add_mat,...
    y,hl, h, T,max_max_l, str_vn_cn, str_cn_vn,dc, dv, theta,nsigma,csigma1,...
    mnk,  cnstl1,gray_labels, code_seq)
q = numel(cnstl1);
M=size(h,1);
N = size(h,2);
[HD0_gf, HD0_cplx] = HD_QAM(y, cnstl1);
Synd0_gf = decod_prod1(HD0_gf',h_arr,M, mul_mat, add_mat);
Synd0_bin = double(Synd0_gf==0);
failed_mx = M-sum(Synd0_bin);

l=0;

seqgf = HD0_gf;

Synd_bin = Synd0_bin;
Synd_bin0 = Synd0_bin;
HD_gf = HD0_gf;

N0 = 2*csigma1^2;
[p_Si_given_y, ~, ~] = Prb_Si_given_y(y, N0, cnstl1);
cmf = p_Si_given_y;
gf_sort = cmf;
b = 1:q;
for i = 1 : N
    [p_Si_given_y(:,i),b] = sort(p_Si_given_y(:,i));
    gf_sort(:,i) = b;
    cmf(:,i) = cumsum(p_Si_given_y(:,i));
end
cc=0;
plt = 0;
while l<T
    if l>-1
        cc=1;
    end
    l = l + 1;

    nse =  cc*nsigma*randn(1,N);


    WSH = Synd_bin*hl;
    E = WSH + nse;
    if plt
        plot(E); ylim([-5 5])
    end
    i2 = find(E<theta);
    if length(i2)>mnk
        [~,i3] = mink(E(i2),mnk);
    else
        i3 = 1:length(i2);
    end
    i1=i2(i3);
    i4 = nan(size(i1));
    for j = 1 :length(i1)
        rnd = rand();
        i4(j) = find(cmf(:,i1(j))>rnd,1);

    end
    for i5 = 1 : length(i1)
        HD_gf(i1(i5)) = gf_sort(i4(i5),i1(i5))-1;
    end

    Synd_gf = decod_prod1(HD_gf',h_arr,M, mul_mat, add_mat);
    Synd_bin = double(Synd_gf==0);
    failed = M-sum(Synd_bin);
    if failed < failed_mx
        Synd_bin0 = Synd_bin;
        seqgf = HD_gf;
        failed_mx = failed;
        if failed==0
            break
        end
    end

end
if failed_mx>0 && max_max_l>l
    fld_idx = find(Synd_bin0==0);
    fld_vns=str_cn_vn(fld_idx);
    combinedArray = [fld_vns{:}]; % Combine all arrays into one
    uniqueIntegers = unique(combinedArray); % Find unique integers
    counts = histc(combinedArray, uniqueIntegers); % Count occurrences of each unique integer
    [~, i1] = sort(counts,'descend');
    vn1 = uniqueIntegers(i1);

    for o=1:min(3, length(vn1))
        dgf = seqgf;
        [~,i2] = sort(abs(y(o)-cnstl1));
        cnddt_gf = i2-1;
        if counts(i1(o))==dv(vn1(o))
            for i = 1:numel(cnddt_gf)
                l=l+1;
                if l > max_max_l
                    break;
                end
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





