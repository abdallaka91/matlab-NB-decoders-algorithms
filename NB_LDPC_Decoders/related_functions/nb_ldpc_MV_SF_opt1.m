function [iter, dec_seq, success_dec] = nb_ldpc_MV_SF_opt1(LLR_2, vi,v0v1, nui,L, iter_max, mul_mat, add_mat, div_mat, combs, h)

v0 = v0v1(1);
v1 = v0v1(2);
M = size(h,1);
N = size(h,2);

q = size(LLR_2,1);
Wmn = cell(M,1);
Dwnm = cell(M,1);
Tml = cell(M,1);
Tn1 = cell(M,1);
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
        L1 = L;%vi^nui;
        mult = -ones(1, dc(i));
        tml = Lnm_3d{i}(:,1);
        coef = h(i,idx1);
        jj = 1:dc(i);
        for ii=1:nui
            jj(jj==idx_chng(ii))=[];
        end
        for j0 = jj
            mult(j0) = mul_mat(tml(j0)+1, coef(j0)+1); %Tml{i}(l,:).*(h_gf(i,idx1));
        end
        mult(idx_chng)=[];
        Ss = sum_arr_gf_dec(mult, add_mat);
        Sl = zeros(L1,1);

        for l = 1 :L1
            Ss1 = Ss;
            tml = Lnm_3d{i}(:,1);
            for j0 = 1 : nui
                II1 = idx_chng(j0);
                II2 = combs(l,j0);
                tml(II1)=Lnm_3d{i}(II1,II2);
                temp = mul_mat(Lnm_3d{i}(II1,II2)+1, coef(II1)+1);
                Ss1 = add_mat(Ss1+1, temp+1);
            end
            Sl(l) = Ss1;
            cc = zeros(1,dc(i));
            for j0 = 1 : dc(i)
                cc(j0) = div_mat(Sl(l)+1, coef(j0)+1);
                cc(j0) = add_mat(cc(j0)+1, tml(j0)+1);
            end
            R{i}(:, l) = cc;
            ii2 = R{i}(:,l)+1;
%             tempp = tml==Lnm_3d{i}(:,1);
            for j1 = 1 : dc(i)
                if l==1 %sum(tempp)-tempp(j1)==dc(i)-1
                    Wmn{i}(j1, ii2(j1)) = Wmn{i}(j1, ii2(j1)) + v0;
                    Wn(idx1(j1), ii2(j1)) = Wn(idx1(j1), ii2(j1)) + v0;
                else
                    Wmn{i}(j1, ii2(j1)) = Wmn{i}(j1, ii2(j1)) + v1;
                    Wn(idx1(j1), ii2(j1)) = Wn(idx1(j1), ii2(j1)) + v1;
                end
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