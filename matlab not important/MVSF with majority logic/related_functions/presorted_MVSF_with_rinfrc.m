function [iters, dec_seq, success_dec] = ...
    presorted_MVSF_with_rinfrc(LLR_2, max_iter, mul_mat, add_mat, div_mat,...
    h,str_cn_vn, dc,str_vn_cn, dv, dev_lsts, nm, v_weights, max_attempt, y_bin_nse)

[M,N] = size(h);
q = size(LLR_2,2);

llr3 = nan(M,q);
gf3=nan(M,q);

relb1 = nan(1,N);

for n=1:N
    [a,b]=sort(LLR_2(n,:), 'descend');
    relb1(n) = a(1)-a(2);
    gf3(n,:)=b-1;
    llr3(n,:)=a;
end
[~,id1] = sort(relb1,'ascend');

gf3_srt = gf3(id1,:);
llr3_srt = llr3(id1,:);

gf3_srt0=gf3_srt;

temp0(id1) = gf3_srt(:,1);
s0 = decod_prod(temp0, h, str_cn_vn, mul_mat,add_mat);
nvld = s0>0;
nvldcnt = sum(nvld);
nvldcnt11=nvldcnt;
for u = 1 : 1
    for n = 1 : N
        temp1 = gf3_srt(:,1);
        for d = 2 :3
            temp1(n) = gf3_srt(n,d);
            temp2(id1) = temp1;
            s1 = decod_prod(temp2, h, str_cn_vn, mul_mat,add_mat);
            nvld1 = s1>0;
            nvldcnt1 = sum(nvld1);
            if nvldcnt1<nvldcnt
                temp3=temp2;
                nvldcnt = nvldcnt1;
                a = gf3_srt(n,1);
                gf3_srt(n,1) = gf3_srt(n,d);
                gf3_srt(n,d) = a;
                llr3_srt(n,1) = llr3_srt(n,1) + 0;
                break
            end
        end
    end
    gf3(id1,:) = gf3_srt;
    llr3(id1,:) = llr3_srt;
    for i = 1 : N
        llr3(i,gf3(i,:)+1) = llr3(i,:);
    end

end

LLR_22=llr3;

[~,HD1] = max(LLR_22,[], 2);
HD1 = HD1'-1;
s2 = decod_prod(HD1, h, str_cn_vn, mul_mat,add_mat);
nvld2 = s2>0;
nvldcnt21 = sum(nvld2);
[iters, dec_seq, success_dec, synd] = presorted_MVSF_4(LLR_2, max_iter, mul_mat, add_mat, div_mat,...
    h,str_cn_vn, dc,str_vn_cn, dv, dev_lsts, nm, v_weights, y_bin_nse);

[iters, dec_seq, success_dec, synd] = presorted_MVSF_4(LLR_22, max_iter, mul_mat, add_mat, div_mat,...
    h,str_cn_vn, dc,str_vn_cn, dv, dev_lsts, nm, v_weights, y_bin_nse);
end
