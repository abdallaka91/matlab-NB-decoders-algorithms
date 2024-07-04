function [v2c_llr, v2c_gf]= truncat_v2c(dc, nm, Mv2c)
v2c_llr = -ones(dc, nm);
v2c_gf = -ones(dc, nm);
    for j = 1 : dc
        [a,b] = sort(Mv2c(j,:));
        v2c_llr(j, :) = a(1:nm);
        v2c_gf(j, :) = b(1:nm) - 1;
    end
end