function [iter, dec_seq, success_dec] = nb_ldpc_MV_SF(LLR_2, vi,v0v1, nui, iter_max, mul_mat, add_mat, div_mat, combs, h)

v0 = v0v1(1);
v1 = v0v1(2);
M = size(h,1);
N = size(h,2);

q = size(LLR_2,1);
Wmn = cell(M,1);
Dwnm = cell(M,1);
Tml = cell(M,1);
R = cell(M,1);
Vmn = cell(M,1);
Qnm_3d = cell(M,1);
Lnm_3d = cell(M,1);
str_cn_vn = cell(M,1);
mn_idx_m = cell(M,1);

dc = zeros(M,1);

Wn = LLR_2';

for i = 1 : M
    str_cn_vn{i, 1} = find(h(i,:));
    dc(i) = length(str_cn_vn{i});
    Tml{i} = zeros(dc(i), size(combs,1));
    R{i} = zeros(dc(i), size(combs,1));
    Vmn{i} = zeros(dc(i),size(combs,1));
end

iter = 0;

for i = 1 : M
    Dwnm{i} = zeros(dc(i),q);
    Wmn{i} = zeros(dc(i),q);
end

while iter<iter_max
    iter = iter+1;
    for i = 1 : M
        idx1 = str_cn_vn{i,1};

        Qnm_3d{i} = zeros(dc(i),q);
        Lnm_3d{i} = zeros(dc(i),vi);
        mn_idx = zeros(dc(i),1);

        for j = 1:dc(i)
            Dwnm{i}(j, :) = Wn(idx1(j),:) - Wmn{i}(j, :);
            temp1 = Dwnm{i}(j, :);
            [a,b] = sort(temp1,'descend');
            %             mn_idx(j0,2) = ;
            mn_idx(j,1) = a(1)-a(2);
            b0 = b - 1;
            Qnm = b0;
            Qnm_3d{i}(j, :) = Qnm;
            Lnm_3d{i}(j, :) = Qnm(1:vi);
            %             Dwnm{i}(j, :) = Dwnm{i}(j, b);%---------------------------------->>>>> hereeee
        end

        [a,b] = sort(mn_idx(:,1));
        idx_chng= b(1:nui);
        mn_idx_m{i} = idx_chng;
        L1 = vi^nui;
        for l = 1 :L1
            temp = Lnm_3d{i};
            for j0 = 1 : nui
                II1 = idx_chng(j0);
                II2 = combs(l,j0);
                temp(II1,1) = temp(II1,II2);
            end
            Tml{i}(:,l) = temp(:,1);
        end

        Sl = zeros(L1,1);
        for l = 1 :L1
            mult = -ones(1, dc(i));
            tml = Tml{i}(:,l);
            coef = h(i,idx1);
            for j0 = 1 : dc(i)
                mult(j0) = mul_mat(tml(j0)+1, coef(j0)+1); %Tml{i}(l,:).*(h_gf(i,idx1));
            end
            Sl(l) = sum_arr_gf_dec(mult, add_mat);
            %             Sl(l)
            cc = zeros(1,dc(i));
            for j0 = 1 : dc(i)
                cc(j0) = div_mat(Sl(l)+1, coef(j0)+1);
                cc(j0) = add_mat(cc(j0)+1, tml(j0)+1);
            end
            R{i}(:, l) = cc;
        end

        %     Vmn{i} is zeros(M,dc(i),size(combs,1));
        for l = 1 :L1
            a1 = Tml{i}(:,l);
            b1 = (squeeze(Qnm_3d{i}(:,1)));
            eq = a1==b1;

            % iiv0 = (eq);
            % iiv1 = (~eq);
            % Vmn{i}(iiv0,l) = v0;
            % Vmn{i}(iiv1,l) = v1;
            if prod(eq)
                Vmn{i}(:,l) = v0;
            else
                Vmn{i}(:,l) = v1;
            end


            ii2 = R{i}(:,l)+1;
            for j1 = 1 : dc(i)
                Wmn{i}(j1, ii2(j1)) = Wmn{i}(j1, ii2(j1)) + Vmn{i}(j1, l);
                Wn(idx1(j1), ii2(j1)) = Wn(idx1(j1), ii2(j1)) + Vmn{i}(j1, l);
            end

        end
    end
    [~, dec_seq] = max(Wn,[],2);
    dec_seq = dec_seq' - 1;

    synd = Inf(1,M);

    for j1 = 1 : M
        idx1 = str_cn_vn{j1};
        tempc = zeros(dc(j1),1);
        for j0 = 1 : dc(j1)
            tempc(j0) = mul_mat(dec_seq(idx1(j0))+1, h(j1,idx1(j0))+1);
        end
        synd(j1) = sum_arr_gf_dec(tempc, add_mat);

    end

    is_not_code = sum(synd)>0;
    success_dec = ~is_not_code;

    if success_dec
        break
    end

end