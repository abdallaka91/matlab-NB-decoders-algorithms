function H = binGaussElim(H)
rows = size(H, 1);
cols = size(H, 2);

r = 1;
for c = cols - rows + 1:cols
    if H(r,c) == 0
        % Swap needed
        for r2 = r + 1:rows
            if H(r2,c) ~= 0
                tmp = H(r, :);
                H(r, :) = H(r2, :);
                H(r2, :) = tmp;
            end
        end

        % Ups...
        if H(r,c) == 0
            error('H is singular');
        end
    end

    % Forward substitute
    for r2 = r + 1:rows
        if H(r2, c) == 1
            H(r2, :) = xor(H(r2, :), H(r, :));
        end
    end

    % Back Substitution
    for r2 = 1:r - 1
        if H(r2, c) == 1
            H(r2, :) = xor(H(r2, :), H(r, :));
        end
    end

    % Next row
    r = r + 1;
end