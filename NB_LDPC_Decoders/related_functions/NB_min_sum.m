function [z, dec_seq, is_code] = NB_min_sum(Niter, LLR_2, dc, lst_comb_all, str_cn_vn, q, h, h_gf)%,max_iter, add_mat, mul_mat, div_mat, v)

is_code = false;
p = log2(q);
M = size(h, 1);
N = size(h,2);
for j1 = 1:q
    h_inf(:,:, j1) = -1e20*h;
end
m_chk_to_var = -1e20*ones(M,N,q);
m_var_to_chk = m_chk_to_var;
m_var_to_chk_t = zeros(M,N,1);
a_posteriori_mat = 0*zeros(q,N);

for z = 1 : Niter
    for jj1 = 1 : q

        if z == 1
            for i = 1 : M
                idx1 = str_cn_vn{i,1};
                aprior = LLR_2(jj1,:);
                m_var_to_chk_t(i, idx1) = aprior(idx1);
            end

        else
            for i = 1 : M
                idx1 = str_cn_vn{i,1};
                aprior = LLR_2(jj1,:);
                m_var_to_chk_t(i, idx1) = aprior(idx1) + a_posteriori_mat(jj1,idx1) -...
                    m_chk_to_var(i,idx1, jj1);
            end

        end
        m_var_to_chk(:,:,jj1) = m_var_to_chk_t;
    end
    m_chk_to_var = h_inf;
    %----------------------------------------

    for i = 1 : M
        idx1 = str_cn_vn{i,1};
        l1 = dc(i)-1;
        lst_comb1 = lst_comb_all{i,1}+1;
        for j2 = 1 : q^l1
            for j1 = 1 : dc(i)
                idx2 = idx1;
                idx2(j1) = [];
                sm = 0;
                idx3 = 1:dc(i);
                idx3(j1) = [];
                idx4 = lst_comb1(j2, idx3);
                for j3 = 1 : l1
                    i2n = idx2(j3);
                    i3n = idx4(j3);
                    sm = sm +  m_var_to_chk(i,i2n, i3n);
                end
                i4n = idx1(j1);
                i5n = lst_comb1(j2, j1);
                if sm>=m_chk_to_var(i,i4n, i5n)
                    m_chk_to_var(i,i4n, i5n)= sm;
                end
            end
        end
    end

    a_posteriori_mat = LLR_2;
    for j1 = 1 : q
        temp = sum(m_chk_to_var(:,:,j1), 1);
        a_posteriori_mat(j1,:) = a_posteriori_mat(j1,:) + temp;
    end

    [~,dec_seq] = max(a_posteriori_mat);
    dec_seq = dec_seq - 1;

    dec_seq_gf = gf(dec_seq,p);

    a = dec_seq_gf*h_gf';
    a1 = a.x;
    if sum(a1)==0
        is_code = 1;
        break
    end
end