function  y = GF_pol_enc_dec_mat_mul(c,G_arr, mul_mat, add_mat)
%     y(h_arr(k)) = add_mat(1+y(h_arr(k)), 1+mul_mat(1+c(h_arr(k+1)), 1+h_arr(k+2)));
y = zeros(size(c));
N3=length(G_arr);
for k = 1 : 3: N3
    b1=mul_mat(c(G_arr(k))+1, G_arr(k+2)+1);
    y(G_arr(k+1))=add_mat(y(G_arr(k+1))+1, b1+1);
end
