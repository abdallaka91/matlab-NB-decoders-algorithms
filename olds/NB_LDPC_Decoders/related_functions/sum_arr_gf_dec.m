function y = sum_arr_gf_dec(x, add_mat)
q = size(add_mat,1);
y=0;
for i = 1 : length(x)
    y = add_mat(y+1, x(i)+1);
end
