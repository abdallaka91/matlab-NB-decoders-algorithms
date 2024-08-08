function [c2v_llr, c2v_gf] = ECN_Fw_Bw_c2v(v2c_gf,v2c_llr, dc,q, nm_v2c, nm_c2v, add_mat, comp_ECN)
c2v_llr = zeros(dc,nm_c2v);
c2v_gf = zeros(dc,nm_c2v);

ECNf_llr = zeros(dc-3,nm_v2c);
ECNf_gf = zeros(dc-3,nm_v2c);

for i = 2 : dc-1
    if i==2
        I(1,:) =v2c_gf(1,:);
        I(2,:) =v2c_llr(1,:);
    else
        I(1,:) = ECNf_gf(i-2,:);
        I(2,:) = ECNf_llr(i-2,:);
    end
    U(1,:) = v2c_gf(i,:);
    U(2,:) = v2c_llr(i,:);
    if i<=dc-2
        V =    ECN4(q, nm_v2c, I, U, add_mat);
        ECNf_gf(i-1,:) = V(1,:);
        ECNf_llr(i-1,:) = V(2,:);
    else
        V =    ECN4_c2v(q, nm_v2c, nm_c2v, I, U, add_mat);
        c2v_llr(dc, :) = V(2,:)-min(V(2,:));
        c2v_gf(dc,:) = V(1,:);
    end
end


ECNb_llr = zeros(dc-3,nm_v2c);
ECNb_gf = zeros(dc-3,nm_v2c);

for i = dc-2 : -1: 1
    if i==dc-2
        I(1,:) =v2c_gf(dc,:);
        I(2,:) =v2c_llr(dc,:);
    else
        I(1,:) = ECNb_gf(i,:);
        I(2,:) = ECNb_llr(i,:);
    end
    U(1,:) = v2c_gf(i+1,:);
    U(2,:) = v2c_llr(i+1,:);
    if i>=2
    V = ECN4(q, nm_v2c, I, U, add_mat);
        ECNb_gf(i-1,:) = V(1,:);
        ECNb_llr(i-1,:) = V(2,:);
    else
        V = ECN4_c2v(q, nm_v2c, nm_c2v, I, U, add_mat);
        c2v_llr(1, :) = V(2,:)-min(V(2,:));
        c2v_gf(1,:) = V(1,:);
    end
end

for i = 1 : dc-2
    if i==1
        I(1,:) =v2c_gf(1,:);
        I(2,:) =v2c_llr(1,:);
    else
        I(1,:) = ECNf_gf(i-1,:);
        I(2,:) = ECNf_llr(i-1,:);
    end
    if i<=dc-3
        U(1,:) = ECNb_gf(i,:);
        U(2,:) = ECNb_llr(i,:);
    else
        U(1,:) = v2c_gf(dc,:);
        U(2,:) = v2c_llr(dc,:);
    end
    V =    ECN4_c2v(q, nm_v2c, nm_c2v, I, U, add_mat);

    c2v_llr(i+1, :) = V(2,:)-min(V(2,:));
    c2v_gf(i+1,:) = V(1,:);

end