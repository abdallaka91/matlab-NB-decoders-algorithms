function inverse_matrix = matrix_inversion_decimal(matrix)
    % Check if the input matrix is square
    [rows, cols] = size(matrix);
    if rows ~= cols
        error('Input matrix must be square for inversion.');
    end
    
    % Augment the matrix with the identity matrix
    augmented_matrix = [matrix, eye(rows)];
    
    % Perform Gauss-Jordan elimination
    for i = 1:rows
        % Partial pivoting to avoid division by zero
        [~, max_row] = max(abs(augmented_matrix(i:rows, i)));
        max_row = max_row + i - 1;
        if max_row ~= i
            augmented_matrix([i max_row], :) = augmented_matrix([max_row i], :);
        end
        
        % Divide the pivot row to make the diagonal element 1
        pivot = augmented_matrix(i, i);
        augmented_matrix(i, :) = augmented_matrix(i, :) / pivot;
        
        % Eliminate non-zero elements below the pivot
        for j = 1:rows
            if j ~= i
                multiplier = augmented_matrix(j, i);
                augmented_matrix(j, :) = augmented_matrix(j, :) - multiplier * augmented_matrix(i, :);
            end
        end
    end
    
    % Extract the inverted matrix from the augmented matrix
    inverse_matrix = augmented_matrix(:, cols+1:end);
end