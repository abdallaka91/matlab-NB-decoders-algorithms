function LLR_gf = LLR_BPSK_GFq_2D(cdata, simga_l)
p = size(cdata,2);
Q = 2^p;
N = size(cdata,1);
LLR_gf = zeros(Q,N);
cdata1 = cdata;% reshape(cdata, p, N)';
data_to_compare = de2bi(0:Q-1, p,2);
for i = 1:N
    input = cdata1(i,:);
    llr_arr = 2*input/(simga_l^2);
    HD = sign(llr_arr);
    HD = -HD/2+0.5;
    conf =abs(llr_arr);
    for j = 1: Q
        XOR  = mod(HD+data_to_compare(j,:), 2);
        llr_o = -sum(XOR.*conf);
        LLR_gf(j,i) = llr_o;
    end
    
end

