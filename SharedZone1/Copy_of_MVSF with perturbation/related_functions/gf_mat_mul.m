function y = gf_mat_mul(a,b,add_mat,mul_mat)

[ma, na] = size(a);
[~, nb] = size(b);

y = zeros(ma, nb);

for i = 1 : ma
    for j = 1 : nb
        rw = a(i,:);
        cl = b(:,j);
        for k=1 : na
            y(i,j) = add_mat(1+y(i,j), 1+mul_mat(1+rw(k),1+cl(k)));
        end
    end
end



