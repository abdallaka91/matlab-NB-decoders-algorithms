clear
% rng(0);
load arith_64;
clc

cnt11=0;
cnt111=0;
EbN0_dB = 0; % SNR per bit in dB
EbN0 = 10^(EbN0_dB/10); % Linear SNR
sigm = sqrt(1/(2*EbN0)); % Noise standard deviation

while cnt11<1000
    nH = 8;
    nL = 20;
    nc = 8;
    nb=5;
    nopM=nL+2;
    cnt11=cnt11+1;

    T = nan(nH, nL);
    TT = nan(nH,nL);
    q = 64;
    p=6;
    words = (0:q-1);
    alph_bin =  fliplr(dec2bin(words, p) - 48);
    alph_bin_mod = (-1).^alph_bin;
    LLR_20 = zeros(1, q);
    y_bin = randi([0 1], 1, 6);
    nse = sigm*randn(size(y_bin));
    y_bin_nse = y_bin + nse;
    a1 = -LLR_simple3(y_bin_nse,1024 , -inf, q,1, alph_bin, LLR_20);
    [a, aa]=sort(round(a1), 'ascend');
    aa = aa-1;

    y_bin_nse = y_bin + nse;
    b1 = -LLR_simple3(y_bin_nse,1024 , -inf, q,1, alph_bin, LLR_20);
    [b, bb]=sort(round(b1), 'ascend');
    bb = bb-1;

    a = a(1:nL);
    a_gf = aa(1:nL);
    b = b(1:nH);
    b_gf = bb(1:nH);
    c=nan(nH, nL);
    cc=c;
    for i = 1 : nH
        for j=1:nL
            cc(i,j)=add_mat(1+aa(i),1+bb(j));
            c(i,j) = a(i)+b(j);
        end
    end

    cnt=0;
    E=[];
    EE=[];%add_mat(aa(1)+1, bb(1)+1);
    T(1:nb,1) = c(1:nb);
    TT(1:nb,1) = cc(1:nb,1);
    srtr = T(1:nb,1);
    srtrg = TT(1:nb,1);
    si = 1:nb;
    sj =ones(1,nb);

    j=1;

    nop=0;


    Ti = false(size(T));
    Ti(1:nb,1)=true;
    H = 1; Hb=~H;
    while nop<nopM
        nop=nop+1;
        [m,n] = min(srtr);
        n1 = srtrg(n);
        i=si(n);
        j=sj(n);
        if sum(EE==srtrg(n))==0
            E=[E m];
            EE=[EE srtrg(n)];
            cnt = cnt+1;
        end
        if i==nH ||j==nL || cnt==nc
            break
        end
        if i==1
            H = 1; Hb=~H;
        end
        if i==nb
            H = 0; Hb=~H;
        end
        if ~Ti(i+Hb, j+H)
            i1 = i+Hb;
            j1 = j+H;
        else
            i1 = i+H;
            j1 = j+Hb;
        end
        Ti(i1, j1) = true;
        T(i1,j1) = a(i1)+b(j1);
        TT(i1,j1) = add_mat(aa(i1)+1, bb(j1)+1);
        srtr(n) = T(i1,j1);
        srtrg(n) = TT(i1,j1);
        si(n) = i1;
        sj(n) = j1;
        if  nop==nopM
            break
        end
    end
    nnn=sort(c(:));
    if length(E)==nc
        ss=E-nnn(1:nc);
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