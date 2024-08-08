function binary_table = deviation_hamm_lst(p, nones)
    binary_table = [];
        for i = 0:(2^p - 1)
        binary_str = dec2bin(i, p);
        num_ones = sum(binary_str == '1');        
        if num_ones <= nones
            binary_word = arrayfun(@(x) str2double(x), binary_str);
            binary_table = [binary_table; binary_word];
        end
    end
end