function [lst2, lst2b] = hmm3(lst,p, df)
lst1=bin2dec(lst);
N=2^p;
lst2=[];
lst2b=[];
% for i=N-1:-1:0
for i=0:N-1
    n=length(lst1);
    dsts=zeros(n,1);
    for j=1:n
        dsts(j)=hamming_distance(i,lst1(j),p);
    end
    if prod(dsts>=df(1) & dsts<=df(2))==1
        lst2=[lst2; i];
        lst2b=[lst2b; dec2bin(i,p)];
    end
end
end
function dist = hamming_distance(x,y, p)
dist=sum(dec2bin(bitxor(x, y),p) == '1');
end