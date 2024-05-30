function [iters,Trial, dec_seq, success_dec] = sweep_corr(LLR_2, mul_mat, add_mat, div_mat, h,str_cn_vn, dc)
M = size(h,1);
N = size(h,2);
q = size(LLR_2,2);
Wmn = cell(M,1);
most_reliabe_idx = zeros(1, N);
second_reliabe_idx = zeros(1, N);
most_reliabe_llr = zeros(1, N);
second_reliabe_llr = zeros(1, N);

for nn = 1 : N
    arr = LLR_2(nn,:);
    [m1,max_idx1] = max(arr);
    most_reliabe_idx(nn) = max_idx1;
    arr(max_idx1) = -inf;
    [m2,max_idx2] = max(arr);
    second_reliabe_idx(nn) = max_idx2;
    most_reliabe_llr(nn) = m1;
    second_reliabe_llr(nn) = m2;
end
relibly = second_reliabe_llr - second_reliabe_llr;
[~, id_sort] = sort(relibly, 'ascend');
h = h(:, id_sort);
LLR_3 = LLR_2(:,id_sort);

dc = zeros(M,1);

for i = 1 : M
    str_cn_vn{i, 1} = find(h(i,:));
    dc(i) = length(str_cn_vn{i});
end