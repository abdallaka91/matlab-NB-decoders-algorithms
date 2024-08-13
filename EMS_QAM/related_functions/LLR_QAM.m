function [LLR_2, HD1] = LLR_QAM(y_cmp, norm_const,sigm, gray_labels)

q = numel(norm_const);
N = numel(y_cmp);
HD1 = zeros(1,N);
LLR_2 = zeros(N,q);

for i = 1 : N
    d = abs(y_cmp(i)-norm_const.');
    [dmn,HD1(i)] = min(d);
    HD1(i) = gray_labels(HD1(i));
    LLR_2(i,:) = (d.^2-dmn^2)*2/sigm^2;
end
end


