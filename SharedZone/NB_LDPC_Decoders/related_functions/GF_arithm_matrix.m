function mat = GF_arithm_matrix(q, arith_func)
p = log2(q);
alph = 0:q-1;
alph_gf = gf(alph, p);
mat = zeros(q,q);
mat_GF = gf(mat, p);

if strcmp(arith_func,'mul')

    for i = 1 : q
        for j = 1 : q
            mat_GF(i, j) = alph_gf(i)*alph_gf(j);
            mat(i, j) = double(mat_GF.x(i, j));
        end
    end
elseif strcmp(arith_func,'add')
    for i = 1 : q
        for j = 1 : q
            mat_GF(i, j) = alph_gf(i)+alph_gf(j);
            mat(i, j) = double(mat_GF.x(i, j));
        end
    end
elseif strcmp(arith_func,'div')
    for i = 1 : q
        for j = 1 : q
            if j~= 1
                mat_GF(i, j) = alph_gf(i)/alph_gf(j);
                mat(i, j) = double(mat_GF.x(i, j));
            else
                mat(i, j) = nan;
            end
        end
    end
else
    mat = [];
end