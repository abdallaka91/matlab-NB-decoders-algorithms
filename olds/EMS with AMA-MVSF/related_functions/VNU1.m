function Mv2c = VNU1(N, dv, Mv2c, Mc2v, APP, str_vn_cn, str_cn_vn)
for j = 1: N
    idx1 = str_vn_cn{j};
    for i = 1 : dv(j)
        i1 = idx1(i);
        idx2 = str_cn_vn{i1};
        i2 = find(idx2==j,1);
        Mv2c{i1}(i2,:) = APP(j,:)- Mc2v{i1}(i2,:);
    end
end