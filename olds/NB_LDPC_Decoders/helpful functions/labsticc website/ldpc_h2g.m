function [H_sys, G] = ldpc_h2g(H, q)


[r, n] = size(H);
k = n - r;

H_gf = gf(H, log2(q)); % 转为GF(q)
H_flr = fliplr(H_gf);
H_flr_rref = rref_GF(H_flr);
H_sys = rot90(H_flr_rref,2);
P = H_sys(:, 1:k);
G = [eye(k), P'];

end
