function decimal_numbers = bi2de(reversed_binary_matrix)

reversed_binary_matrix = fliplr(reversed_binary_matrix);

% Convert each row of the binary matrix to a string
binary_strings = num2str(reversed_binary_matrix);

% Convert binary strings to decimal numbers
decimal_numbers = bin2dec(binary_strings);