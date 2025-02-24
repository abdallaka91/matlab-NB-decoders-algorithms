clear
hamming_distance = @(a, b) arrayfun(@(x) sum(dec2bin(bitxor(a, x)) == '1'), b);

N=16;
q=N;
p=log2(N);

t2db=gray_code(p);
t1d=bi2de(t2db);

t2d=t1d;
t2d(1:N/2)=t1d(1:2:N);
t2d(N/2+1:3*N/4-1)=t1d(4:2:N/2);
t2d(1+3*N/4:N-1)=flip(t1d(N/2+4:2:N));
t2d(3*N/4)=t1d(2);
t2d(N)=t1d(N/2+2);


cmbs = nchoosek(1:q, 2);
lc=size(cmbs,1);

hmdst1=zeros(lc,1);
for i=1:lc
    hmdst1(i)=hamming_distance(t1d(cmbs(i,1)), t1d(cmbs(i,2)));
end

hmdst2=zeros(lc,1);
for i=1:lc
    hmdst2(i)=hamming_distance(t2d(cmbs(i,1)), t2d(cmbs(i,2)));
end

dst1=1;

vv=[hmdst1 hmdst2];
idx1=hmdst1==1;
s1=sum(hmdst2(idx1));
idx2=hmdst1==2;
s2=sum(hmdst2(idx2));
vv1=[hmdst1(idx1) hmdst2(idx1)];
vv2=[hmdst1(idx2) hmdst2(idx2)];