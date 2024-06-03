function [iters, dec_seq, success_dec, syndrms1] = ...
    presorted_MVSF_with_rinfrc1(LLR_2, max_iter, mul_mat, add_mat, div_mat,...
    h,str_cn_vn, dc, dev_lsts, nm, v_weights, max_attempt)

[M,N] = size(h);
iters = 0;
LLR_20 = LLR_2;
tt =0;
syndrms1 = nan(M,max_attempt*max_iter);

while tt<max_attempt
    id1 = zeros(M,1);
    s1 = zeros(1,N);
    tt=tt+1;
    [iters0, dec_seq, success_dec, syndrms] = presorted_MVSF_5(LLR_2, max_iter, mul_mat, add_mat, div_mat,...
        h,str_cn_vn, dc, dev_lsts, nm, v_weights);
    % syndrms1(:,(tt-1)*max_iter+1:tt*max_iter)=syndrms;
    synd = syndrms(:,end);
    iters = iters+iters0;

    if success_dec
        break
    else
        for i = 1 : M
            if ~synd(i)
                aa = syndrms(i,:);
                i1 = find(aa==1,1,'last');
                if isempty(i1)
                    i1=0;
                end
                if i1 < max_iter/2
                    id1(i)=2;
                else
                    id1(i)=1;
                end

            end
        end
        LLR_2 = LLR_20;
        for i = 1 : M
            s1(str_cn_vn{i}) = id1(i);
        end
        for j = 1 : N
            if s1(j)==2
                LLR_2(j,dec_seq(j)+1) = 1;
                LLR_2(j,:) = LLR_2(j,:) -1;
            % elseif s1(j)==1
            %     LLR_2(j,dec_seq(j)+1) = 1;
            %     LLR_2(j,:) = LLR_2(j,:) - 1;
            end
        end
    end
end