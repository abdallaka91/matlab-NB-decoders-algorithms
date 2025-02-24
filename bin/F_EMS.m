function output_LLR = F_EMS(input_LLR1, input_LLR2, h, table, GF, nm, bubbleregion, offset, l)
    output_LLR = ones(1, GF) * 1e3;
    
    sorted_GF1 = insertion_sort(input_LLR1, GF);
    sorted_GF2 = insertion_sort(input_LLR2, GF);
    
    if input_LLR1(sorted_GF1(3)) > input_LLR2(sorted_GF2(3))
        for i = 1:nm
            for j = 1:nm
                if bubbleregion(i, j) == 1
                    temp = input_LLR1(sorted_GF1(i)) + input_LLR2(sorted_GF2(j));
                    res_sym = table.ADDGF(sorted_GF1(i), sorted_GF2(j));
                    output_LLR(res_sym) = min(output_LLR(res_sym), temp);
                end
            end
        end
    else
        for i = 1:nm
            for j = 1:nm
                if bubbleregion(i, j) == 1
                    temp = input_LLR1(sorted_GF1(j)) + input_LLR2(sorted_GF2(i));
                    res_sym = table.ADDGF(sorted_GF1(j)+1, sorted_GF2(i)+1);
                    output_LLR(res_sym) = min(output_LLR(res_sym), temp);
                end
            end
        end
    end
    
    best_GF = insertion_sort(output_LLR, GF)-1;
    tmp_LLR = ones(1, GF) * (output_LLR(best_GF(nm)+1) + offset);
    
    for q = 1:nm
        tmp_LLR(best_GF(q)+1) = output_LLR(best_GF(q)+1);
    end
    
    output_LLR = tmp_LLR;
end

function [sorted_indices, cnt] = insertion_sort(LLRs, array_size)
    sorted_indices = zeros(1, array_size);
    cnt = 1;
    for i = 1:array_size
        sorted_indices(i) = i;
    end
    for i = 2:array_size
        j = i;
        while (j > 1) && (LLRs(sorted_indices(j - 1)) > LLRs(i))
            sorted_indices(j) = sorted_indices(j - 1);
            j = j - 1;
        end
        sorted_indices(j) = i;
        
    end

end
