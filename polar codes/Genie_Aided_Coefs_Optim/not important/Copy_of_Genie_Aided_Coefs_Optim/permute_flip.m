function [y, hmdst2, t2d]=permute_flip(bits,cmbs,idx, t2d_bin, p1, q)
bits=round(bits);
for i=1:q
    if bits(i)==1
        t2d_bin(i,:)=-t2d_bin(i,:)+1;
    end
end
t2d=bi2de(t2d_bin);

lc=length(idx);
hmdst2=zeros(lc,1);
for i=1:lc
    hmdst2(i)=hamming_distance(t2d(cmbs(idx(i),1)), t2d(cmbs(idx(i),2)));
end

y=sum(hmdst2==p1);

end

function y=hamming_distance(a,b)
y= sum(dec2bin(bitxor(a, b)) == '1');
end