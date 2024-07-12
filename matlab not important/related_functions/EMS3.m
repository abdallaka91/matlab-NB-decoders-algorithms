function [iters, dec_seq, success_dec, synd,APP] = EMS3(LLR_2, max_iter, mul_mat, add_mat, div_mat,...
    h,str_cn_vn, dc,str_vn_cn,dv, nm_v2c, nm_c2v, c2v_comp_fact, comp_ECN)
% c2v_comp = 50.33;
M = size(h,1);
N = size(h,2);
q = size(LLR_2,2);
success_dec = 0;
iters = 0;

Mv2c = cell(M,1);
Mc2v = cell(M,1);

for i = 1:M
    Mv2c{i} = zeros(dc(i),q);
    Mc2v{i} = zeros(dc(i),q);
end

iter = 0;
APP = LLR_2;
while iter<max_iter
    iter = iter+1;
    Mv2c = VNU1(N, dv, Mv2c, Mc2v, APP, str_vn_cn, str_cn_vn);
    APP = LLR_2;
    iters = iters+1;
    for i = 1 : M
        temp1 = Mv2c{i};
        [v2c_llr, v2c_gf]= truncat_v2c(dc(i), nm_v2c, temp1);
        v2c_gfp = v2c_gf;

        idx1 = str_cn_vn{i};
        coefs = h(i, idx1);
        for j = 1 : dc(i)
            for k = 1 : nm_v2c
                v2c_gfp(j,k) = mul_mat(v2c_gf(j, k)+1, coefs(j)+1);
            end
        end

        [c2v_llr, c2v_gfp] = ECN_Fw_Bw_c2v(v2c_gfp,v2c_llr, dc(i),q, nm_v2c, nm_c2v, add_mat,comp_ECN);
        %          [c2v_llr11, c2v_gfp11] = ECN_Fw_Bw_1111(v2c_gfp,v2c_llr, dc(i),nm, add_mat);
        %         for j = 1 : dc(i)
        %             for k = 1 : nm_v2c
        %                 if isnan(c2v_gfp(j, k))
        %                     for iii = 0:q-1
        %                         cnd1 = true;
        %                         for iiii=1:nm_c2v
        %                             if v2c_gf(j, k)==iii
        %                                 cnd1 = false;
        %                                 break
        %                             end
        %                         end
        %                         if cnd1
        %                             break
        %                         end
        %                     end
        %                     c2v_gfp(j, k) =  iii;
        %                     c2v_llr(j, k) = max(c2v_llr(j,:)) +c2v_comp_fact;
        %                 end
        %             end
        %         end

        for j = 1 : dc(i)
            for kk=nm_v2c:nm_c2v
                if isnan(c2v_llr(j,kk))
                    break
                end
            end
            temp = c2v_llr(j,kk)*ones(1,q)+c2v_comp_fact;
            for k = 1 : kk
                temp_c2v_gf = div_mat(c2v_gfp(j, k)+1, coefs(j)+1);
                temp(temp_c2v_gf+1) = c2v_llr(j,k);
            end
            Mc2v{i}(j,:) = temp;
            APP(idx1(j),:) = APP(idx1(j),:) + temp;
        end

    end
    %     APP = APP_U(N, dv, Mc2v, LLR_2, str_vn_cn, str_cn_vn);
    [mm, dec_seq] = min(APP,[],2);
    APP = bsxfun(@minus, APP, mm);
    dec_seq = dec_seq' - 1;
    sum(dec_seq~=0);
    synd=decod_prod(dec_seq, h, str_cn_vn, mul_mat, add_mat);
    is_not_code = sum(synd)>0;
    success_dec = ~is_not_code;
    if success_dec
        break
    end

end

