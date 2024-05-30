function [iters,Trial, dec_seq, success_dec, Wn] = presorted_MVSF_3(LLR_2, max_iter, mul_mat, add_mat, div_mat,...
    h,str_cn_vn, dc, dev_lsts,dev_pos, nm, v_weights, Dev_pos_cnt, PERM_rng, PERM, max_trial)
M = size(h,1);
N = size(h,2);
q = size(LLR_2,2);
p=log2(q);
Wmn = cell(M,1);
most_reliabe_idx = zeros(1, N);
second_reliabe_idx = zeros(1, N);
most_reliabe_llr = zeros(1, N);
second_reliabe_llr = zeros(1, N);

LLR_3 = LLR_2;

for nn = 1 : N
    arr = LLR_3(nn,:);
    [tt1, max_idx1] = max(arr);
    arr(max_idx1) = -inf;
    most_reliabe_llr(nn) = tt1;
    most_reliabe_idx(nn) = max_idx1;
    Q = [1:max_idx1-1 max_idx1+1:q];
    whgt_fact = Q;
    i0=0;
    for i = Q
        i0=i0+1;
        a = de2bi(max_idx1-1, p);
        b = de2bi(i-1, p);
        n1 = sum(a~=b);
        if n1==1
            whgt_fact(i0)=1;
        elseif n1==2
            whgt_fact(i0) = 1.0;
        else
            whgt_fact(i0)=1.0;
        end
    end
    arr(Q) = arr(Q).*whgt_fact;
    [tt2, max_idx2] = max(arr);
    a = de2bi(max_idx1-1, p);
    b = de2bi(max_idx2-1, p);
    n1 = sum(a~=b);
    second_reliabe_llr(nn) = tt2;
    second_reliabe_idx(nn) = max_idx2;
end

relibly = most_reliabe_llr - second_reliabe_llr;
[~, id_sort] = sort(relibly, 'ascend');
success_dec = 0;
iters = 0;
Wn = LLR_2;
for pp = 1 : max_trial
    for i = 1 : M
        Wmn{i} = zeros(dc(i),q);
    end
    Mv2c_LLR_T = cell(M,1);
    Mv2c_GF_T = cell(M,1);

    for i = 1:M
        Mv2c_GF_T{i} = zeros(dc(i),nm);
        Mv2c_LLR_T{i} = 0;
    end

    iter = 0;
    while iter<max_iter
        iter = iter+1;
        iters = iters+1;
        for i = 1 : M
            idx1 = str_cn_vn{i,1};
            for j = 1:dc(i)
                Dwnm = Wn(idx1(j),:) - Wmn{i}(j, :);
                [a,b] = sort(Dwnm, 'descend');
                Mv2c_LLR_T{i}(j,:) = a(1)-a(2);
                Mv2c_GF_T{i}(j,:) = b(1:nm)-1;
            end

            idx1 = str_cn_vn{i,1};
            Mv2c_LLR_i1 = Mv2c_LLR_T{i};
            Mv2c_GF_i1 = Mv2c_GF_T{i};
            U1 = Mv2c_LLR_i1;
            [~,idxx2] = sort(U1, 'descend');
            Mv2c_GF_i = Mv2c_GF_i1(idxx2,:);
            coef1 = h(i, idx1);
            coef = coef1(idxx2);
            L = size(dev_lsts{i},1);
            syndrom_i_GF = zeros(L,1);
            lst_comb_gf = zeros(L, dc(i));
            dev_lsts_i = dev_lsts{i};
            % dev_pos_i = dev_pos{i};
            dev_pos_i = dev_lsts_i-1;
            temp_gf = zeros(1,dc(i));

            for l = 1 : L
                dev_i_l = dev_lsts_i(l,:);
                for j1 = 1 : dc(i)
                    temp_gf(j1) = Mv2c_GF_i(j1, dev_i_l(j1));
                    lst_comb_gf(l,j1) = temp_gf(j1);
                    temp_gf(j1) = mul_mat(temp_gf(j1)+1,coef(j1)+1);
                end
                syndrom_i_GF(l) = sum_arr_gf_dec(temp_gf, add_mat);
            end

            nm_c2v = L;
            msg_gf_i2j0 = nan(dc(i),nm_c2v);
            msg_gf_i2j = msg_gf_i2j0;
            msg_llr_i2j0 = nan(dc(i),nm_c2v);
            msg_llr_i2j = msg_llr_i2j0;
            for l = 1 : L
                temp_lst_comb_gf1 = lst_comb_gf(l,:);
                for addr = 1 : dc(i)
                    temp_gf1 = div_mat(syndrom_i_GF(l)+1, coef(addr)+1);
                    temp_gf1 = add_mat(temp_gf1+1, temp_lst_comb_gf1(addr)+1);
                    msg_gf_i2j0(addr,l) = temp_gf1;
                    msg_llr_i2j0(addr,l) = v_weights(2);
                    if l==1
                        msg_llr_i2j0(addr,1)=v_weights(1);
                    end

                end
                msg_gf_i2j(idxx2,:) = msg_gf_i2j0;
                msg_llr_i2j(idxx2,:) = msg_llr_i2j0;
            end
            for j1 = 1 : dc(i)
                for l = 1 : nm_c2v
                    Wmn{i}(j1, msg_gf_i2j(j1,l)+1) = Wmn{i}(j1, msg_gf_i2j(j1,l)+1) + msg_llr_i2j(j1,l);
                    Wn(idx1(j1), msg_gf_i2j(j1,l)+1) = Wn(idx1(j1),msg_gf_i2j(j1,l)+1) + msg_llr_i2j(j1,l);
                end
            end
        end

        [mm, dec_seq] = max(Wn,[],2);
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
    % if ~success_dec
    %     max_iter = 25;
    %     Wn = LLR_2;
    %     rndperm = randperm(PERM_rng,PERM);
    %     most_reliabe_llr1 = most_reliabe_llr(id_sort);
    %     second_reliabe_llr1 = second_reliabe_llr(id_sort);
    %     temp = most_reliabe_llr1(rndperm);
    %     most_reliabe_llr1(rndperm) = second_reliabe_llr1(rndperm);
    %     second_reliabe_llr1(rndperm) = temp;
    %     for uu = 1 : PERM
    %         Wn(id_sort(uu),second_reliabe_idx(id_sort(uu))) = second_reliabe_llr1(uu);
    %         Wn(id_sort(uu), most_reliabe_idx(id_sort(uu))) = most_reliabe_llr1(uu);
    %     end
    % else
    %     break;
    % end
end
Trial = pp;
