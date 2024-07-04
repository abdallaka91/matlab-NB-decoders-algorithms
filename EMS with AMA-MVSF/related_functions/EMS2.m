function [iters, dec_seq, success_dec, synd,APP] = EMS2(LLR_2, max_iter, mul_mat, add_mat, div_mat,...
    h,str_cn_vn, dc,str_vn_cn,dv, nm, nc, c2v_comp_fact)
% c2v_comp = 50.33;
M = size(h,1);
N = size(h,2);
q = size(LLR_2,2);
success_dec = 0;
iters = 0;

Mv2c = cell(M,1);
Mv2c_gf_perm = cell(M,1);
Mc2v = cell(M,1);
c2v_gf = cell(M,1);

v2c_llr_cell = cell(M,1);
v2c_gf_cell = cell(M,1);

for i = 1:M
    v2c_llr_cell{i} = zeros(dc(i), nm);
    v2c_gf_cell{i} = v2c_llr_cell{i};
    Mv2c{i} = zeros(dc(i),q);
    Mv2c_gf_perm{i} = zeros(dc(i),q);
    Mc2v{i} = zeros(dc(i),q);
    c2v_gf{i} = zeros(dc(i),nm);
end

for i = 1 : M
    idx1 = str_cn_vn{i};
    for j = 1 : dc(i)
        for k = 1 : q
            Mv2c_gf_perm{i}(j,k) = mul_mat(k, 1+h(i, idx1(j)));
        end
    end
end


iter = 0;


while iter<max_iter
    iter = iter+1;
    Mv2c = VNU1(N, dv, Mv2c, Mc2v, LLR_2, str_vn_cn, str_cn_vn);
    iters = iters+1;
    for i = 1 : M
        [v2c_llr, v2c_gfp]= truncat_v2c1(dc(i), nm, Mv2c{i}, Mv2c_gf_perm{i}, v2c_llr_cell{i}, v2c_gf_cell{i});
        [c2v_llr, c2v_gfp] = ECN_Fw_Bw(v2c_gfp,v2c_llr, dc(i), nm, nc, add_mat);
        for j = 1 : dc(i)
            for k = 1 : nm
                c2v_gf{i}(j,k) = div_mat(c2v_gfp(j, k)+1, h(i, str_cn_vn{i}(j))+1);
            end
        end
        for j = 1 : dc(i)
            Mc2v{i}(j,:) = repmat(max(c2v_llr(j,:))+c2v_comp_fact,1,q);
            for k = 1 : nm
                Mc2v{i}(j,c2v_gf{i}(j,k)+1) = c2v_llr(j,k);
            end

        end
    end
    APP = APP_U(N, dv, Mc2v, LLR_2, str_vn_cn, str_cn_vn);

    [~, dec_seq] = min(APP,[],2);
    mm=min(APP,[],2);
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

