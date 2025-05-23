function [iters, dec_seq, success_dec, synd,Wn] = presorted_MVSF_4(LLR_2, max_iter, mul_mat, add_mat, div_mat,...
    h,str_cn_vn, dc,str_vn_cn, dv, dev_lsts, nm, v_weights, y_bin_nse)
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



% function [iters, dec_seq, success_dec, synd,Wn] = presorted_MVSF_4(LLR_2, max_iter, mul_mat, add_mat, div_mat,...
%     h,str_cn_vn, dc,str_vn_cn, dv, dev_lsts, nm, v_weights, y_bin_nse)
% 
% bb=40;
% M = size(h,1);
% N = size(h,2);
% q = size(LLR_2,2);
% p=log2(q);
% words = (0:q-1);
% alph_bin =  fliplr(dec2bin(words, p) - 48);
% alph_bin_mod = (-1).^alph_bin;
% 
% Wmn = cell(M,1);
% Dwnm = cell(M,1);
% 
% success_dec = 0;
% iters = 0;
% Wn = LLR_2;
% 
% for i = 1 : M
%     Wmn{i} = zeros(dc(i),q);
%     Dwnm{i} = zeros(dc(i),q);
% end
% Mv2c_LLR_T = cell(M,1);
% Mv2c_GF_T = cell(M,1);
% msg_c2v_gf = cell(M,1);
% msg_c2v_llr = cell(M,1);
% 
% 
% for i = 1:M
%     msg_c2v_gf{i} = nan(dc(i),size(dev_lsts{i},1));
%     msg_c2v_llr{i} = nan(dc(i),size(dev_lsts{i},1));
%     Mv2c_GF_T{i} = zeros(dc(i),nm);
%     Mv2c_LLR_T{i} = 0;
% end
% 
% iter = 0;
% while iter<max_iter
%     iter = iter+1;
%     iters = iters+1;
%     for i = 1 : M
%         idx1 = str_cn_vn{i,1};
%         for j = 1:dc(i)
%             Dwnm{i}(j, :) = Wn(idx1(j),:) - Wmn{i}(j, :);
%             [a,b] = sort(Dwnm{i}(j, :), 'descend');
%             Mv2c_LLR_T{i}(j,:) = a(1)-a(2);
%             Mv2c_GF_T{i}(j,:) = b(1:nm)-1;
%         end
% 
%         idx1 = str_cn_vn{i,1};
%         Mv2c_LLR_i1 = Mv2c_LLR_T{i};
%         Mv2c_GF_i1 = Mv2c_GF_T{i};
%         U1 = Mv2c_LLR_i1;
%         [~,idxx2] = sort(U1, 'descend');
%         Mv2c_GF_i = Mv2c_GF_i1(idxx2,:);
%         coef1 = h(i, idx1);
%         coef = coef1(idxx2);
%         L = size(dev_lsts{i},1);
%         syndrom_i_GF = zeros(L,1);
%         lst_comb_gf = zeros(L, dc(i));
%         dev_lsts_i = dev_lsts{i};
%         temp_gf = zeros(1,dc(i));
% 
%         for l = 1 : L
%             dev_i_l = dev_lsts_i(l,:);
%             for j1 = 1 : dc(i)
%                 temp_gf(j1) = Mv2c_GF_i(j1, dev_i_l(j1));
%                 lst_comb_gf(l,j1) = temp_gf(j1);
%                 temp_gf(j1) = mul_mat(temp_gf(j1)+1,coef(j1)+1);
%             end
%             syndrom_i_GF(l) = sum_arr_gf_dec(temp_gf, add_mat);
%         end
% 
%         nm_c2v = L;
%         msg_gf_i2j0 = nan(dc(i),1);
% 
%         for l = 1 : L
%             temp_lst_comb_gf1 = lst_comb_gf(l,:);
%             for addr = 1 : dc(i)
%                 temp_gf1 = div_mat(syndrom_i_GF(l)+1, coef(addr)+1);
%                 temp_gf1 = add_mat(temp_gf1+1, temp_lst_comb_gf1(addr)+1);
%                 msg_gf_i2j0(addr) = temp_gf1;
%             end
%             msg_c2v_gf{i}(idxx2,l) = msg_gf_i2j0;
%         end
% 
%         idx1 = str_cn_vn{i,1};
%         for j1 = 1 : dc(i)
%             i11=1;
%             for l = 1 : nm_c2v
%                 if l>1
%                     i11=2;
%                 end
%                 Wmn{i}(j1, msg_c2v_gf{i}(j1,l)+1) = Wmn{i}(j1, msg_c2v_gf{i}(j1,l)+1) + v_weights(i11);
%                 Wn(idx1(j1), msg_c2v_gf{i}(j1,l)+1) = Wn(idx1(j1),msg_c2v_gf{i}(j1,l)+1) + v_weights(i11);
%             end
%         end
%     end
% 
%     [mm, dec_seq] = max(Wn,[],2);
%     Wn = bsxfun(@minus, Wn, mm);
%     dec_seq = dec_seq' - 1;
%     sum(dec_seq~=0);
%     synd=decod_prod(dec_seq, h, str_cn_vn, mul_mat, add_mat);
%     is_not_code = sum(synd)>0;
%     success_dec = ~is_not_code;
%     if success_dec
%         break
%     end
% 
%     h11 = zeros(size(h));
% 
%     for i = 1 :M
%         idx1 = str_cn_vn{i,1};
%         if synd(i)==0
%             h11(i,idx1) = 1;
%         else
%             h11(i,idx1) = 0;
%         end
%     end
% 
%      
% end
% 
