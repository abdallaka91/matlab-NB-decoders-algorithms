function [G,Hbar] = Generator_matrix_G_from_full_rank_H(H, add_mat, mul_mat, div_mat)
[M,N]=size(H);
K=N-M;
n_k = M;
H1 = H(:,1:K);
H2 = H(:,1+K:end);

H2_inv = matrix_inversion_GF(H2, mul_mat, add_mat, div_mat);

Hbar = [gf_mat_mul(H2_inv, H1,add_mat,mul_mat) eye(n_k)];
G = [eye(K) (gf_mat_mul(H2_inv, H1,add_mat,mul_mat))'];
end
function inverse_matrix = matrix_inversion_GF(matrix, mul_mat, add_mat, div_mat)
[rows, cols] = size(matrix);

augmented_matrix = [matrix, eye(rows)];

for i = 1:rows
    [~, max_row] = max(abs(augmented_matrix(i:rows, i)));
    max_row = max_row + i - 1;
    if max_row ~= i
        augmented_matrix([i max_row], :) = augmented_matrix([max_row i], :);
    end
    pivot = augmented_matrix(i, i);
    augmented_matrix(i, :) = div_mat(1+augmented_matrix(i, :), 1+pivot);

    for j = 1:rows
        if j ~= i
            multiplier = augmented_matrix(j, i);
            ml1 = mul_mat(1+multiplier,1+augmented_matrix(i, :));
            for k = 1 : length(augmented_matrix)
                augmented_matrix(j, k) = add_mat(1+augmented_matrix(j, k),1+ml1(k));
            end
        end
    end
end

inverse_matrix = augmented_matrix(:, cols+1:end);
end

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
end