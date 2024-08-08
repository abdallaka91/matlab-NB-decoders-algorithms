function [iters0, dec_seq0, success_dec0] = ...
    presorted_MVSF_with_rinfrc1(LLR_1, max_iter, mul_mat, add_mat, div_mat,...
    h,str_cn_vn, dc,str_vn_cn, dv, dev_lsts, nm, v_weights, max_attempt, y_bin_nse)

[iters0, dec_seq0, success_dec0, synd0,LLR_2] = presorted_MVSF_4(LLR_1, max_iter, mul_mat, add_mat, div_mat,...
    h,str_cn_vn, dc,str_vn_cn, dv, dev_lsts, nm, v_weights, y_bin_nse);
% tr1 = 1;
% Ad = 10;
% Rng = 60;
% success_dec0=0;
% if ~success_dec0 && tr1
%     LLR_2 = LLR_1;
%     [~,HD0] = max(LLR_2,[], 2);
%     HD0 = HD0'-1;
%     HD1 = HD0;
%     s1 = decod_prod(HD0, h, str_cn_vn, mul_mat,add_mat);
%     nvld1 = s1>0;
%     nvldcnt2 = sum(nvld1);
%     nmx=nvldcnt2;
% 
%     [M,N] = size(h);
%     q = size(LLR_2,2);
% 
% 
%     for d = 1 :10
%         nmx=nvldcnt2;
%         LLR_22=LLR_2;
%         nflp=0;
% 
%         relbl = zeros(1,N);
%         for n=1:N
%             [a,b] = maxk(LLR_22(n,:),2);
%             relbl(n)=a(1)-a(2);
%         end
%         [relbl_srt, id1] = sort(relbl);
% 
% 
%         for n=randperm(40,10)
%             [a,b] = maxk(LLR_22(id1(n),:),2);
%             LLR_22(id1(n),b(1))=a(2)+0.1;
%             %                 LLR_22(id1(n),b(2))=a(1);
%             HD1(id1(n)) = b(2)-1;
%             s2 = decod_prod(HD1, h, str_cn_vn, mul_mat,add_mat);
%             nvld2 = s2>0;
%             nvldcnt21 = sum(nvld2);
% 
%                             if nvldcnt21<=nmx-2
%                                 nflp=nflp+1;
%                                 nmx = nvldcnt21;
%                             else
%                                 HD1(id1(n)) = b(1)-1;
%                                 LLR_22(id1(n),b(1))=a(1);
%                                 LLR_22(id1(n),b(2))=a(2);
%                             end
%         end
%         [~,HD1] = max(LLR_22,[], 2);
%         HD1 = HD1'-1;
%         s2 = decod_prod(HD1, h, str_cn_vn, mul_mat,add_mat);
%         nvld2 = s2>0;
%         nvldcnt21 = sum(nvld2);
% 
%         aa=LLR_2(id1,:);
%         bb=LLR_22(id1,:);
% 
% 
%         [iters0, dec_seq0, success_dec0, synd0] = presorted_MVSF_4(LLR_22, max_iter, mul_mat, add_mat, div_mat,...
%             h,str_cn_vn, dc,str_vn_cn, dv, dev_lsts, nm, v_weights, y_bin_nse);
%         if success_dec0
%             break
%         end
% 
%     end
% end

