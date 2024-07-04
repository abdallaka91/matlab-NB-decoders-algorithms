function [v2c_llr, v2c_gfp]= truncat_v2c1(dc, nm, Mv2c,Mv2c_gf_perm,  v2c_llr, v2c_gfp)
    for j = 1 : dc
        [v2c_llr(j, :),b] = mink(Mv2c(j,:), nm);
        v2c_gfp(j, :) = Mv2c_gf_perm(j,b);
    end
end