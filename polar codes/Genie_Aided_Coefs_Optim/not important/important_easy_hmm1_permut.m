clear

% This code efficiently permutes the sequence elements from 0 to 2^p - 1,
% where p is a strictly even integer.
% The permutation ensures that every pair in the original sequence with a
% Hamming distance of 1 (or p-1) is mapped to a new pair with a 
% Hamming distance of p-1 (or 1), respectively.

% z1 = [0 1 1 0];
% Z=[0 1 1 0];
% O=[1 0 0 1];



hamming_distance = @(a, b) arrayfun(@(x) sum(dec2bin(bitxor(a, x)) == '1'), b);
p=8;
q=2^p;
t1d=0:q-1;

cmbs = nchoosek(1:q, 2);
lc=size(cmbs,1);

hmdst1=zeros(lc,1);
for i=1:lc
    hmdst1(i)=hamming_distance(t1d(cmbs(i,1)), t1d(cmbs(i,2)));
end
%%
z1 = [0 1 1 0];
Z=[0 1 1 0];
O=[1 0 0 1];

for i=1:log(q)/log(4)-1
z1 = kron(z1, O) + kron(1 - z1, Z);  
end

t1d_bin=fliplr(de2bi(t1d,p));
t2d_bin=t1d_bin;

for i=1:q
    if z1(i)==1
        t2d_bin(i,:)=-t2d_bin(i,:)+1;
    end
end

t2d=bi2de(t2d_bin);

cmbs = nchoosek(1:q, 2);
lc=size(cmbs,1);

hmdst2=zeros(lc,1);
for i=1:lc
    hmdst2(i)=hamming_distance(t2d(cmbs(i,1)), t2d(cmbs(i,2)));
end

vv=[hmdst1 hmdst2];
idx1=hmdst1==1;
s1=sum(hmdst2(idx1));
idx2=hmdst1==2;
s2=sum(hmdst2(idx2));
vv1=[hmdst1(idx1) hmdst2(idx1)];
%%
idx1=find(hmdst1==1);
p1=p-1;
z11=z1;
z11(4)=1;
y=permute_flip(z11,cmbs,idx1, t1d_bin, p1,q);

% vv1=[hmdst1(idx1) hmdst2(idx1)];

