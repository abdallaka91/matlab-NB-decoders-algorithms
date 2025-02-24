clear
hamming_distance = @(a, b) arrayfun(@(x) sum(dec2bin(bitxor(a, x)) == '1'), b);
p=6;
nc=2;
p1=p-nc;
[~, aa] = hmm3("000000" ,p, [p1 p1]);
[~, bb] = hmm3("000001" ,p, [p1 p1]);
aa=string(aa);
bb=string(bb);
cmbs=nchoosek(1:length(aa),2);
dsts=nan(size(cmbs,1),1);
ddd=["000000"];
aa=[aa;bb];
for i=1:5%length(aa)
    cmbs1=nchoosek(1:length(aa), i);
    for j=1:size(cmbs1,1)
        tmp=aa(cmbs1(j,:));
        [~, dd] = hmm3( tmp,p, [p1 p1]);
        dd=string(dd);
        if length(dd)==2
            uu=sum(strcmp(dd(2),aa));
            if ~uu
                aa=[aa; dd(2)];
            end
        end
    end
end


