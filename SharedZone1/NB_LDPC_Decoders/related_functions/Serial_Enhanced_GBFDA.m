function [iter, cn, is_code] = Serial_Enhanced_GBFDA(LLR_2, dc, lst1, q, H,max_iter, ...
    add_mat, mul_mat, div_mat, v)



M = size(H, 1);
N = size(H, 2);
Wmn = cell(M,1);
syndrm = zeros(1,M);
Wn = LLR_2;

for i = 1 : M
    Wmn{i} = zeros(q,dc(i));
end

for iter = 1 : max_iter
    for m = 1 : M
        idx1 = lst1{m};
        Qn = zeros(1,dc(m));
        for j = 1:dc(m)
            temp1 = Wn(:,idx1(j)) - Wmn{m}(:,j);
            [~,rec_info_seq_bit] = max(temp1);
            Qn(j) = rec_info_seq_bit-1;
        end

        s = 0;
        for i = 1 : dc(m)
            temp =  mul_mat(Qn(i)+1, H(m, idx1(i))+1);
            s = add_mat(s+1, temp+1);
        end
        Rn = Qn;
        for i = 1 : dc(m)
            temp = div_mat(s+1, H(m, idx1(i))+1);
            Rn(i) = add_mat(temp+1, Qn(i)+1);
        end

        i = 1;
        for n = idx1
            Wmn{m}(Rn(i)+1,i) = Wmn{m}(Rn(i)+1,i) + v;
            Wn(Rn(i)+1,n) = Wn(Rn(i)+1,n)  + v;
            i = i+1;
        end
    end
    [~,cn] = max(Wn, [],1); cn = cn-1;

    for m = 1 : M
        idx1 = lst1{m};
        s = 0;
        for i = 1 : dc(m)
            temp =  mul_mat(cn(idx1(i))+1, H(m, idx1(i))+1);
            s = add_mat(s+1, temp+1);
        end
        syndrm(m) = s;
    end

    is_code = sum(syndrm)==0;
    if is_code
        break
    end
end