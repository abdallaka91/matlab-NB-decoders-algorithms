function [K, Z, success_dec] = P_SFDP_func(y_bin_nse, alphb_bin, iter_max, mul_mat, add_mat, div_mat, h, list_CN, list_VN, dc, etta)


alphb = bi2de(alphb_bin);
M = size(h,1);
N = size(h,2);
q = size(alphb_bin,1);
p = log2(q);


Z = nan(1, N);
Z_bin = ones(size(y_bin_nse));
Z_bin(y_bin_nse<0) = 0;

for i = 1 : N
Z(i) = bi2de( Z_bin(i,:));
end

K = 0;

EXIij = zeros(size(h));

Ej = zeros(N, q-1);
Vj = zeros(N, q-1);
Vjmax = zeros(N, 1);
jk1 = 1:N;
ik = 1:M;
pp = -1;
scnd = false;
while K<iter_max
    K = K+1;

    for i = ik
        idx1 = list_CN{i,1};
        Qn = Z(idx1);
        s = 0;
        for m = 1 : dc(i)
            temp =  mul_mat(Qn(m)+1, h(i, idx1(m))+1);
            s = add_mat(s+1, temp+1);
        end
        for m = 1 : dc(i)
            temp = div_mat(s+1, h(i, idx1(m))+1);
            EXIij(i, idx1(m)) = add_mat(temp+1, Qn(m)+1);
        end
    end
    
    
    for j = jk1
        idx1 = list_VN{j};
        Zj = Z(j);
        Zj_bin = Z_bin(j,:);
        Zj_cand = [0:Zj-1 Zj+1:q-1];
        Zj_cand_bin = alphb_bin(Zj_cand+1,:);
        Vj(j,:) = Zj_cand;
        bao_zj = bina_asym_op(Zj_bin, y_bin_nse(j,:));
        sigmj = EXIij(idx1, j);
        nZj = sum(Zj==sigmj);
        Ej(j, :) = - bao_zj - etta(nZj+1);

        for i = 1 : q-1
            Zj_c = Zj_cand(i);
            Zj_c_bin = Zj_cand_bin(i,:);
            bao_zj_c = bina_asym_op(Zj_c_bin, y_bin_nse(j,:));
            nZj_c = sum(Zj_c==EXIij(idx1, j));
            Ej(j, i) = Ej(j, i) + bao_zj_c + etta(nZj_c+1);               ;
        end
    end
    [Ejmax, i1] = max(Ej,[],2);
    

    for i2 = 1 : N
        Vjmax(i2) = Vj(i2,i1(i2));
    end

    [~, jjj] = maxk(Ejmax, 2);
    jjj1 = jjj(1);
    jjj2 = jjj(2);

    if (jjj1)~=pp
        jk = jjj1;
        pp = jjj1;
    else
        scnd = true;
        jk = jjj2;
        pp = jjj2;
    end
        
    Z(jk) = Vjmax(jk);
    Z_bin(jk,:) = de2bi(Vjmax(jk), p);
    ik = list_VN{jk}';
    
    jk1 = zeros(1,sum(dc(ik)));

    tt = 1;
    for j1 = 1 :  length(ik)
        jk1(tt:tt+dc(j1)-1) = list_CN{ik(j1)};
        tt = tt + dc(j1);
    end
%     plot(Z)
    
    synd = Inf(1,M);

    for j1 = 1 : M
        idx1 = list_CN{j1};
        tempc = zeros(dc(j1),1);
        for j0 = 1 : dc(j1)
            tempc(j0) = mul_mat(Z(idx1(j0))+1, h(j1,idx1(j0))+1);
        end
        synd(j1) = sum_arr_gf_dec(tempc, add_mat);

    end

    is_not_code = sum(synd)>0;
    success_dec = ~is_not_code;

    if success_dec
        break
    end

end
end

function x = bina_asym_op(z, y)
x = sum((2*z-1).*y);
end

