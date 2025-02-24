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