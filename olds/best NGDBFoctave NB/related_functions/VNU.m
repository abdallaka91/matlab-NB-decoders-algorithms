function Mv2c = VNU(N, dv, Mv2c, Mc2v, APP, str_vn_cn, str_cn_vn)
for j = 1: N
    idx1 = str_vn_cn{j};
    for i = 1 : dv(j)
        i1 = idx1(i);
        idx2 = str_cn_vn{i1};
        i2 = find(idx2==j,1);
        Mv2c{i1}(i2,:) = APP(j,:)- Mc2v{i1}(i2,:);
    end
end

% function Mv2c = VNU(N, dv, Mv2c, Mc2v, LLR_2, str_vn_cn, str_cn_vn)
% for j = 1: N
%     idx1 = str_vn_cn{j};
%     for i = 1 : dv(j)
%         i1 = idx1(i);
%         idx2 = str_cn_vn{i1};
%         i2 = find(idx2==j,1);
%         Mv2c{i1}(i2,:) = LLR_2(j,:);
%         idx11 = idx1([1:i-1 i+1:dv(j)]);
%         for i3 = 1  : dv(j)-1
%             i31 = idx11(i3);
% 
%             idx21 = str_cn_vn{i31};
%             i22 = find(idx21==j,1);
%             Mv2c{i1}(i2,:) = Mv2c{i1}(i2,:) + Mc2v{i31}(i22,:);
%         end
% 
%     end
% end