function [iter, dec_seq, success_dec] = presorted_MVSF(LLR_2, max_iter, mul_mat, add_mat, div_mat, h,str_cn_vn, dc, ...
    dev_lsts,dev_pos, nm, Ofst)
% Ofst = 0.7;
M = size(h,1);
N = size(h,2);
q = size(LLR_2,2);

Wn = LLR_2;
% Wn = Wn/max(max(Wn));
Wmn = cell(M,1);
for i = 1 : M
    Wmn{i} = zeros(dc(i),q);
end


Mv2c_LLR = cell(M,1);
Mv2c_LLR_T = cell(M,1);
Mv2c_GF_T = cell(M,1);

for i = 1:M
    Mv2c_LLR{i} = zeros(dc(i),q);
    Mv2c_GF_T{i} = zeros(dc(i),nm);
    Mv2c_LLR_T{i} = zeros(dc(i),nm);
end

iter = 0;

mmx = 1;
while iter<max_iter
    iter = iter+1;
    for i = 1 : M
        idx1 = str_cn_vn{i,1};
        for j = 1:dc(i)
            Dwnm = Wn(idx1(j),:) - Wmn{i}(j, :)/mmx;
            [a,b] = sort(Dwnm);
            Mv2c_LLR_T{i}(j,:) = a(1:nm)-a(1);
            Mv2c_GF_T{i}(j,:) = b(1:nm)-1;
        end
    end


    for i = 1 : M
        idx1 = str_cn_vn{i,1};
        L = length(dev_lsts{i});
        syndrom_i_LLR = zeros(L,1);
        syndrom_i_GF = zeros(L,1);
        Mv2c_LLR_i = Mv2c_LLR_T{i};
        Mv2c_GF_i = Mv2c_GF_T{i};
        lst_comb_gf = zeros(L, dc(i));
        lst_comb_llr = zeros(L, dc(i));
        dev_lsts_i = dev_lsts{i};
        dev_pos_i = dev_pos{i};
        coef = h(i, idx1);
        for l = 1 : L
            temp_llr = zeros(1,dc(i));
            temp_gf = zeros(1,dc(i));
            dev_i_l = dev_lsts_i(l,:);
            for j1 = 1 : dc(i)
                temp_llr(j1) = Mv2c_LLR_i(j1,dev_i_l(j1));
                temp_gf(j1) = Mv2c_GF_i(j1, dev_i_l(j1));
                temp_gf(j1) = mul_mat(temp_gf(j1)+1,coef(j1)+1);
            end
            syndrom_i_LLR(l) = sum(temp_llr);
            syndrom_i_GF(l) = sum_arr_gf_dec(temp_gf, add_mat);
            lst_comb_gf(l,:) = temp_gf;
            lst_comb_llr(l,:) = temp_llr;
        end
        [syndrom_i_LLR1, idxx] = sort(syndrom_i_LLR);
%         syndrom_i_GF1 = syndrom_i_GF(idxx);
%         lst_comb_gf1 = lst_comb_gf(idxx, :);
%         dev_pos_i1 = dev_pos_i(idxx,:);
        for j1 = 1 : dc(i)
            nm1 = 1;
            msg_gf_i2j = zeros(1,nm);
            msg_llr_i2j = zeros(1,nm);
            for l = 1 : L
                dev_pos_i11 = dev_pos_i(idxx(l),:);
                if dev_pos_i11(j1)==0
                    temp_lst_comb_gf1 = lst_comb_gf(idxx(l), :);
                    temp_gf1 = add_mat(syndrom_i_GF(idxx(l))+1, temp_lst_comb_gf1(j1)+1);
                    temp_gf1 = div_mat(temp_gf1+1, coef(j1)+1);
                    if l == 1
                        msg_gf_i2j(nm1) = temp_gf1;
                        msg_llr_i2j(nm1) = syndrom_i_LLR1(l);
                        nm1=nm1+1;
                    elseif sum(temp_gf1==msg_gf_i2j(1:nm1-1))==0
                        msg_gf_i2j(nm1) = temp_gf1;
                        msg_llr_i2j(nm1) = syndrom_i_LLR1(l);
                        nm1=nm1+1;
                    end
                    if nm1>nm
                        break
                    end
                end
            end
            a = 0:q-1;
            least_reliable_GF = setdiff(a, msg_gf_i2j);
            kij = msg_llr_i2j(end) + Ofst;
            Wmn{i}(j1, msg_gf_i2j+1) = msg_llr_i2j;
            Wn(idx1(j1), msg_gf_i2j+1) = Wn(idx1(j1),msg_gf_i2j+1) + msg_llr_i2j;
            if length(least_reliable_GF)>=1
                Wmn{i}(j1, least_reliable_GF+1) = kij;
                Wn(idx1(j1), least_reliable_GF+1) = Wn(idx1(j1),least_reliable_GF+1) + kij;
            end
        end
    end
    [mm, dec_seq] = min(Wn,[],2);
    dec_seq = dec_seq' - 1;
    Wn = bsxfun(@minus, Wn, mm);
%     mmx = max(max(Wn));
%     Wn = Wn/mmx;

    synd=decod_prod(dec_seq, h, str_cn_vn, mul_mat, add_mat);

    is_not_code = sum(synd)>0;
    success_dec = ~is_not_code;

    if success_dec
        break
    end

end