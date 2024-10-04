function [iters, dec_seq, success_dec, synd,APP] = EMS2(LLR_2, max_iter, mul_mat, add_mat, div_mat,...
    h,str_cn_vn, dc,str_vn_cn,dv, nm, nc, c2v_comp_fact, comp_ECN, gray_labels)
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
        [v2c_llr, v2c_gf]= truncat_v2c(dc(i), nm, temp1);
        v2c_gfp = v2c_gf;
        c2v_gf = v2c_gf;
        idx1 = str_cn_vn{i};
        coefs = h(i, idx1);
        for j = 1 : dc(i)
            for k = 1 : nm
                v2c_gfp(j,k) = mul_mat(v2c_gf(j, k)+1, coefs(j)+1);
            end
        end

        [c2v_llr, c2v_gfp] = ECN_Fw_Bw(v2c_gfp,v2c_llr, dc(i),q, nm, nc, add_mat,comp_ECN);
%          [c2v_llr11, c2v_gfp11] = ECN_Fw_Bw_1111(v2c_gfp,v2c_llr, dc(i),nm, add_mat);
        for j = 1 : dc(i)
            for k = 1 : nm
                c2v_gf(j,k) = div_mat(c2v_gfp(j, k)+1, coefs(j)+1);
            end
        end
        for j = 1 : dc(i)
            temp = repmat(max(c2v_llr(j,:))+c2v_comp_fact,1,q);
            for k = 1 : nm
                temp(c2v_gf(j,k)+1) = c2v_llr(j,k);
            end
            Mc2v{i}(j,:) = temp;
            APP(idx1(j),:) = APP(idx1(j),:) + temp;
        end
    end
%     APP = APP_U(N, dv, Mc2v, LLR_2, str_vn_cn, str_cn_vn);
    [mm, dec_seq] = min(APP,[],2);
    dec_seq = gray_labels(dec_seq);
    APP = bsxfun(@minus, APP, mm);
    dec_seq = dec_seq;
    sum(dec_seq~=0);
    synd=decod_prod(dec_seq, h, str_cn_vn, mul_mat, add_mat);
    is_not_code = sum(synd)>0;
    success_dec = ~is_not_code;
    if success_dec
        break
    end

end

