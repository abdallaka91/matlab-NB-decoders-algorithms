function lst_cndtt = Lst_cndtt(y, dev_hamm,pw)

HD1 = sign(y);
cnddts = repmat(HD1, size(dev_hamm,1), 1);
lst_cndtt = nan(1,size(dev_hamm,1));
for ii = 1 : size(dev_hamm,1)
    cnddts(ii,dev_hamm(ii,:))=- HD1(dev_hamm(ii,:));
    lst_cndtt(ii) = 0.5*(1-cnddts(ii,:))*pw;
end
corr = cnddts*y';

[corr_srt, i11] = sort(corr,'descend');
lst_cndtt = lst_cndtt(i11);