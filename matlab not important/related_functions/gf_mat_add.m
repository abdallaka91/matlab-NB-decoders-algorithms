function y = gf_mat_add(a,b,add_mat)

[ma, na] = size(a);

y = zeros(ma, na);

for i = 1 : ma
    for j = 1 : na
        y(i,j) = add_mat(a(i,j)+1, b(i,j)+1);
    end
end