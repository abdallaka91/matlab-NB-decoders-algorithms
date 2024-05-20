clear
N = 1024;
K=4;
Ec = 1;
N0 = 0.5;
Z = exp(-Ec/N0);

L = log2(N);

nn=[0 1];
nnn=nn;
for i = 1 : L-1
    nnn = combvec(nn,nnn);
end
nnn=fliplr(nnn');
mmm=nnn;
for i = 1 : N
    if nnn(i,1)==1
        mmm(i,1)=2*Z-Z^2;
    else
        mmm(i,1)=Z^2;
    end
end
for j = 2 : L
    for i = 1 : N
        if nnn(i,j)==1
            mmm(i,j)=2*mmm(i,j-1)-mmm(i,j-1)^2;
        else
            mmm(i,j) = mmm(i,j-1)^2;
        end
    end
end
Leaves = mmm(:,L);
[Leaves1,I]=sort(Leaves); 
