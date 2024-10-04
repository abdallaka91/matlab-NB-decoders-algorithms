function LLR_gf = LLR_simple4(y, p, simga_l, data_to_compare, alph_bin_mod)
N = size(y,1);
Q = 2^p;
LLR_gf = zeros(Q,N);
cdata1 = y;% reshape(cdata, p, N)';
for i = 1:N
    input = cdata1(i,:);
    llr_arr = 2*input/(simga_l^2);
    HD = sign(llr_arr);
    HD = -HD/2+0.5;
    conf =abs(llr_arr);
    for j = 1: Q
        XOR  = HD~=data_to_compare(j,:);
        llr_o = sum(XOR.*conf);% minus ti make them all negative with max =0 for the HD symbol
        LLR_gf(j,i) = llr_o;
    end
end

