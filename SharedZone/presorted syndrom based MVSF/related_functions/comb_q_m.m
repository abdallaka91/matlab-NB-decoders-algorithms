function y = comb_q_m(current_arr, q, m)

    current_arr(1) = current_arr(1) + 1;

    for j = 1:m
        if current_arr(j)==q
            current_arr(j) = 0;
            current_arr(j+1)=current_arr(j+1)+1;
        else
            break
        end
    end
    y = current_arr;
