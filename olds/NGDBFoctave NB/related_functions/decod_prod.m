function  y = decod_prod(c,h,str_cn_vn,dc, mul_mat, add_mat)
M = size(h,1);
y = zeros(1,M);
% for i = 1 : M
%     idx1 = str_cn_vn{i};
%     vni = c(idx1);
%     coef = h(i, idx1);
%     tempp = 0;
%     for j = 1 : length(idx1)
%         tempp = add_mat(tempp+1, (mul_mat(vni(j)+1, coef(j)+1))+1);
%     end
%     y(i) = tempp;
% end
for i = 1 : M
    tempp = 0;
    for j = 1 : dc(i)
        tempp = add_mat(tempp+1,1+mul_mat(1+c(str_cn_vn{i}(j)), 1+h(i, str_cn_vn{i}(j))));
    end
    y(i) = tempp;
end
