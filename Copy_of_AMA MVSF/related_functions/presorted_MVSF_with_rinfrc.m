function [iters, dec_seq, success_dec, FER_iter, BER_iter, gen_iter] = ...
    presorted_MVSF_with_rinfrc(LLR_2,hd0, max_iter, mul_mat, add_mat, div_mat,...
    h,str_cn_vn, dc, dev_lsts, nm, v_weights, max_attempt,code_seq, FER_iter, BER_iter,gen_iter)

[M,N] = size(h);
iters = 0;
LLR_20 = LLR_2;
tt =0;
iter00 = 0;
while tt<max_attempt
    s1 = zeros(1,N);
    tt=tt+1;
    [iters0, dec_seq, success_dec, synd,~, FER_iter, BER_iter, gen_iter] = presorted_MVSF_4(LLR_2, max_iter, mul_mat, add_mat, div_mat,...
    h,str_cn_vn, dc, dev_lsts, nm, v_weights,code_seq, FER_iter, BER_iter,gen_iter,  iter00);
    iter00 = iter00 + max_iter;
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
%                 b = LLR_2(j,:);
%                 [~,i1] = max(b);
                bb=LLR_2(j,dec_seq(j)+1);
                LLR_2(j,dec_seq(j)+1) = 1;
%                 LLR_2(j,hd0(j)) = 2*bb;
            end
        end
    end
end