function APP = APP_U(N, dv, Mc2v, LLR_2, str_vn_cn, str_cn_vn)
APP = 0*LLR_2;
for j = 1: N
    idx1 = str_vn_cn{j};
    for i = 1 : dv(j)
        APP(j,:) = LLR_2(j,:);
        idx11 = idx1;
        for i3 = 1  : dv(j)-1
            i31 = idx11(i3);
            idx21 = str_cn_vn{i31};
            i22 = find(idx21==j,1);
            APP(j,:) = APP(j,:) + Mc2v{i31}(i22,:);
        end

    end
end