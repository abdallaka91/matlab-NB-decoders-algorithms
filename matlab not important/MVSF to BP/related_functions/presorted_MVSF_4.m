function [iters, dec_seq, success_dec, synd,Wn] = presorted_MVSF_4(LLR_2, max_iter, mul_mat, add_mat, div_mat,...
    h,str_cn_vn, dc, dev_lsts, dv, str_vn_cn,nm, v_weights, cmp_c2v, nc2v)
M = size(h,1);
N = size(h,2);
q = size(LLR_2,2);
p=log2(q);
Wmn = cell(M,1);

success_dec = 0;
iters = 0;
Wn = LLR_2;

for i = 1 : M
    Wmn{i} = zeros(dc(i),q);
end
Mv2c_LLR_T = cell(M,1);
Mv2c_GF_T = cell(M,1);

for i = 1:M
    Mv2c_GF_T{i} = zeros(dc(i),nm);
    Mv2c_LLR_T{i} = zeros(dc(i),nm);
end

iter = 0;
while iter<max_iter

    iter = iter+1;
    iters = iters+1;
    for i = 1 : M
        idx1 = str_cn_vn{i,1};
        for j = 1:dc(i)
            Dwnm = Wn(idx1(j),:) - Wmn{i}(j, :);
            [a,b] = mink(Dwnm, nm);
            a = a - a(1);
            Mv2c_LLR_T{i}(j,:) = a(1:nm);
            Mv2c_GF_T{i}(j,:) = b(1:nm)-1;
        end
        U1 = Mv2c_LLR_T{i}(:,2)-Mv2c_LLR_T{i}(:,1);
        [~,idxx2] = sort(U1, 'descend');
        Mv2c_GF_i = Mv2c_GF_T{i}(idxx2,:);
        Mv2c_llr_i = Mv2c_LLR_T{i}(idxx2,:);
        coef  = h(i, idx1(idxx2));
        L = size(dev_lsts{i},1);
        syndrom_i_GF = zeros(L,1);
        syndrom_i_llr = zeros(L,1);
        lst_comb_gf = zeros(L, dc(i));
        lst_comb_llr = lst_comb_gf;
        temp_gf = zeros(1,dc(i));

        for l = 1 : L
            for j1 = 1 : dc(i)
                temp_gf(j1) = Mv2c_GF_i(j1, dev_lsts{i}(l,j1));
                lst_comb_gf(l,j1) = temp_gf(j1);
                lst_comb_llr(l,j1) = Mv2c_llr_i(j1, dev_lsts{i}(l,j1));
                temp_gf(j1) = mul_mat(temp_gf(j1)+1,coef(j1)+1);
            end
            syndrom_i_GF(l) = sum_arr_gf_dec(temp_gf, add_mat);
            syndrom_i_llr(l) = sum(lst_comb_llr(l,:));
        end

        [syndrom_i_llr, idx3] = sort(syndrom_i_llr,'ascend');
        syndrom_i_GF = syndrom_i_GF(idx3);
        lst_comb_llr = lst_comb_llr(idx3, :);
        lst_comb_gf = lst_comb_gf(idx3, :);

        msg_gf_i2j = nan(dc(i),l);
        for l = 1 : L
            temp_lst_comb_gf1 = lst_comb_gf(l,:);
            for addr = 1 : dc(i)
                temp_gf1 = div_mat(syndrom_i_GF(l)+1, coef(addr)+1);
                temp_gf1 = add_mat(temp_gf1+1, temp_lst_comb_gf1(addr)+1);
                msg_gf_i2j(idxx2(addr),l) = temp_gf1;

            end
        end
        lst_comb_llr(:,idxx2) = lst_comb_llr;
        llr_to_add = nan(dc(i), q);
        for j1 = 1 : dc(i)
            kk=0;
            for l = 1 : L
                if lst_comb_llr(l,j1)==0 && isnan(llr_to_add(j1, msg_gf_i2j(j1,l)+1))
                    llr_to_add(j1, msg_gf_i2j(j1,l)+1)=syndrom_i_llr(l);
                    kk=kk+1;
                    if kk>=nc2v
                        break;
                    end
                end
            end

            mn1 = max(llr_to_add(j1,~isnan(llr_to_add(j1,:))));
            llr_to_add(j1,isnan(llr_to_add(j1,:))) = round(mn1)+cmp_c2v;
            llr_to_add(j1,:)=llr_to_add(j1,:)-min(llr_to_add(j1,:));
        end

        for j1 = 1 : dc(i)
            Wmn{i}(j1, :) =  llr_to_add(j1,:);
            Wn(idx1(j1),:) = Wn(idx1(j1),:) + llr_to_add(j1,:);
        end
    end

    [mm, dec_seq] = min(Wn,[],2);
    Wn = bsxfun(@minus, Wn, mm);
    dec_seq = dec_seq' - 1;
    sum(dec_seq~=0);
    synd=decod_prod(dec_seq, h, str_cn_vn, mul_mat, add_mat);
    is_not_code = sum(synd)>0;
    success_dec = ~is_not_code;
    if success_dec
        break
    end

end

