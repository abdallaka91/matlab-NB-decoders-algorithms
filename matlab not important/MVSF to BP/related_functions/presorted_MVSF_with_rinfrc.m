function [iters, dec_seq, success_dec] = ...
    presorted_MVSF_with_rinfrc(LLR_2, max_iter, mul_mat, add_mat, div_mat,...
    h,str_cn_vn, dc, dev_lsts,dv, str_vn_cn, nm, v_weights, max_attempt, cmp_c2v, nc2v)

[M,N] = size(h);
iters = 0;
LLR_20 = LLR_2;
tt =0;

while tt<max_attempt
    s1 = zeros(1,N);
    tt=tt+1;
    [iters0, dec_seq, success_dec, synd] = presorted_MVSF_4(LLR_2, max_iter, mul_mat, add_mat, div_mat,...
    h,str_cn_vn, dc, dev_lsts, dv, str_vn_cn,nm, v_weights, cmp_c2v, nc2v);
    iters = iters+iters0;
    if success_dec
        break
    else
        LLR_2 = LLR_20;
        for i = 1 : M
            if synd(i)
                s1(str_cn_vn{i}) = 1;
            end
        end
        for j = 1 : N
            if ~s1(j)
                b = LLR_2(j,:);
                [~,i1] = min(b);
                bb=LLR_2(j,dec_seq(j)+1);
                LLR_2(j,dec_seq(j)+1) = 0;
                LLR_2(j,i1) = 2*bb;
            end
        end
    end
end