clear
% rng(0);
% load RNG.mat
% rng(s);
load arith_64;
clc

cnt11=0;
cnt111=0;
nH = 8;
nL = 8;
nc = 8;
nb=8;
zC = 3;
nopM=8*2;
p=6;
q = 2^p;
words = (0:q-1);
alph_bin =  fliplr(dec2bin(words, p) - 48);
alph_bin_mod = (-1).^alph_bin;

sigm = 0.7;
mm=[];
kk=[];
while cnt11<1000
    cnt11=cnt11+1;

    y_bin1 = randi([0 1], 1, p);
    y_bin1_m = 1-2*y_bin1;
    a_dec = bi2de(y_bin1);
    nse = sigm*randn(size(y_bin1_m));
    y_bin_nse = y_bin1_m + nse;
    LLR_20 = zeros(1, q);
    a1 = -LLR_simple3(y_bin_nse,sigm, alph_bin);
    [a, aa]=sort(a1, 'ascend');
    a_gf = aa-1;

    y_bin2 = randi([0 1], 1, p);
    y_bin2_m = 1-2*y_bin2;
    b_dec = bi2de(y_bin2);
    nse = sigm*randn(size(y_bin2_m));
    y_bin_nse = y_bin2_m + nse;
    LLR_20 = zeros(1, q);
    b1 = -LLR_simple3(y_bin_nse,sigm, alph_bin);
    [b, bb]=sort(b1, 'ascend');
    b_gf = bb-1;
    c=nan(nH, nL);
    c_gf=c;

    for i = 1 : nH
        for j=1:nL
            c_gf(i,j)=add_mat(1+a_gf(i),1+b_gf(j));
            c(i,j) = a(i)+b(j);
        end
    end

    [c1,i11]=sort(c(:));
%     C = [c_gf(i11) c1];
    kkk=0;
    C=[c_gf(i11(1)) 0];
    for ii=2:length(c1)
        if sum(c_gf(i11(ii))==C(:,2))==0
            kkk=kkk+1;
            C = [C; [c_gf(i11(ii)) c1(ii)]];
        end
    end

    E1= Lbubble2(a, b, a_gf, b_gf, add_mat, nL, nH, nopM, nb,nc, zC, C);


    
    E=E1(:,2)';
    nnn=C(:,2);
    if length(E)==nc
        ss=E-nnn(1:nc)';
        if max(abs(ss))>1e-5
            
            [d1,d2] = find(ss>1e-5,1);
            mm=[mm d2];
            ff = nc;
            ss=E(1:ff)-nnn(1:ff)';
            if max(abs(ss))>1e-5
                fprintf("a:   " +num2str(E) +'\n');
                fprintf("b:   " +[num2str(nnn(1:nc)')]+'\n--------------------------------------------------------------------------------------\n');
            end
            cnt111=cnt111+1;
        end
    else
        cnt111=cnt111+1;
        kk=[kk;length(E)];
    end
end