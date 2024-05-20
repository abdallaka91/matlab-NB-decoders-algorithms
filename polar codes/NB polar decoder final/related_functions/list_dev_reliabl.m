function [lst_deviation_lst, lst_dev_pos,deviation_lst1,dev_pos1] = list_dev_reliabl(dc_lst, di_lst) % y least 3 di [d0 d1 d2]

dcc = sum(dc_lst);
zones_cnt = length(dc_lst);
lst_deviation_lst = cell(zones_cnt,1);
lst_dev_pos = lst_deviation_lst;

for z = 1 : zones_cnt
    di = di_lst{z};
    dc = dc_lst(z);
    l1 = length(di);
    % D_sz = 1 + dc*di(2) + nchoosek(dc, di(3))*di(3)^2;
    card_Di = zeros(1, l1);
    for i = 0 : l1-1
        if i == 0
            card_Di(i+1) = 1;
        elseif di(i+1)>=1
            card_Di(i+1) = nchoosek(dc, i)*(di(i+1)^i);
        else
            card_Di(i+1) = 0;
        end
    end
    D_sz = sum(card_Di);
    deviation_lst = ones(D_sz, dc);
    dev_pos = zeros(D_sz, dc);
    l = 2;
    for i = 1 : length(di)-1

        comb0 = 2:di(i+1)+1;
        comb00 = comb0;
        for h = 1:i-1
            comb00 = CombVec(comb00,comb0);
        end
        comb00 = comb00';
        cmb = nchoosek(1:dc, i);

        for j = 1 : size(cmb,1)
            idx = ones(1, dc);
            for k = 1 : size(comb00,1)
                idx(cmb(j,:)) = comb00(k,:);
                dev_pos(l, cmb(j,:)) = 1;
                deviation_lst(l,:) = idx;
                l=l+1;
            end
        end
    end

    lst_deviation_lst{z} = deviation_lst;
    lst_dev_pos{z} = dev_pos;
end

combs_cnt = size(lst_deviation_lst{1},1);
comb_zones_lst = 1:combs_cnt;
for i = 2 : zones_cnt
    comb_zones_lst = CombVec(comb_zones_lst,1:size(lst_deviation_lst{i},1));
    combs_cnt = combs_cnt * size(lst_deviation_lst{i},1);
end
comb_zones_lst = comb_zones_lst';

deviation_lst1 = zeros(combs_cnt,dcc);
dev_pos1 = zeros(combs_cnt,dcc);
for i = 1 : combs_cnt
    temp0 = comb_zones_lst(i,:);
    temp1 = [];
    temp2 = [];
    for j = 1 : zones_cnt
        temp1 = [temp1 lst_deviation_lst{j}(temp0(j),:)];
        temp2 = [temp2 lst_dev_pos{j}(temp0(j),:)];
    end
    deviation_lst1(i,:) = temp1;
    dev_pos1(i,:) = temp2;
end
end

