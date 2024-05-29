function [deviation_lst, dev_pos] = syndrome_based_list_dev(dc, di) % y least 3 di [d0 d1 d2]

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
        comb00 = combvec(comb00,comb0);
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