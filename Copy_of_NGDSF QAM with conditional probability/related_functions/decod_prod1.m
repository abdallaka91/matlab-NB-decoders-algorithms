function  y = decod_prod1(c,h_arr,M, mul_mat, add_mat)
y = zeros(1,M);
O=length(h_arr);
for k = 1 : 3: O
    y(h_arr(k)) = add_mat(1+y(h_arr(k)), 1+mul_mat(1+c(h_arr(k+1)), 1+h_arr(k+2)));
end
