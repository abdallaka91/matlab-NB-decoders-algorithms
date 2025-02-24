clear
% rng(0);
load arith_64;
clc

cnt11=0;
cnt111=0;
nH = 8;
nL = 20;
nc = 20;
nb=4;
nopM=nL+2;

while cnt11<1000
    cnt11=cnt11+1;

    T = nan(nH, nL);
    TT = nan(nH,nL);

    a = randperm(50, nH);
    b = randperm(50, nL);
    a=a-min(a);
    b=b-min(b);

    a=sort(a);
    b=sort(b);


    a_gf = randperm( 64, nH)-1;
    b_gf = randperm(64, nL)-1;

    [E, E_gf] = Lbubble(a,b,a_gf, b_gf, add_mat, nL, nH, nopM, nb,nc);

    a=sort(a);
    b=sort(b);
    c=nan(nH, nL);
    c_gf=c;
    for i = 1 : nH
        for j=1:nL
            c_gf(i,j)=add_mat(1+a_gf(i),1+b_gf(j));
            c(i,j) = a(i)+b(j);
        end
    end

    nnn=sort(c(:));
    if length(E)==nc
        ss=E-nnn(1:nc)';
        if max(abs(ss))>1e-5
            E;
            fprintf("a:   " +num2str(E) +'\n');
            fprintf("b:   " +[num2str(nnn(1:nc)')]+'\n--------------------------------------------------------------------------------------\n');

            cnt111=cnt111+1;
        end
    else
        cnt111=cnt111+1;
    end
end