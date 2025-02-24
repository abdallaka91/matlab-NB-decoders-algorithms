clear
% rng(0);
load arith_64;

nm = 8;
nc = 8;
nb=4;

T = nan(nm);
TT = nan(nm);

% a = randperm(50, nm);
% b = randperm(50, nm);
% a = [0 6 13 17 21];
% b = [0 7 15 21 25];
% a=[0    15    29    38    40 42 44 48];
% b=[0    13    30    36    37 39 41 45];

a = randn(1, nm);
b = randn(1, nm);

aa = randperm( 64, nm)-1;
bb = randperm(64, nm)-1;


a=a-min(a);
b=b-min(b);

a=sort(a);
b=sort(b);
c=nan(nm);
cc=c;
for i = 1 : nm
    for j=1:nm
        cc(i,j)=add_mat(1+aa(i),1+bb(j));
        c(i,j) = a(i)+b(j);
    end
end

cnt=0;
E=[];
EE=[];%add_mat(aa(1)+1, bb(1)+1);
T(1:nb,1) = a(1:nb);
TT(1:nb,1) = c(1:nb,1);
srtr = T(1:nb,1);
srtrg = TT(1:nb,1);
si = 1:nb;
sj =ones(1,nb);

j=1;

nop=0;
nopM=100;

Ti = false(size(T));
Ti(1:nb,1)=true;
H = 1; Hb=~H;
while nop<nopM
    nop=nop+1;
    [m,n] = min(srtr);
    n1 = srtrg(n);
    i=si(n);
    j=sj(n);
    E=[E m];
    EE=[EE srtrg(n)];
    cnt = cnt+1;
    if i==nm ||j==nm || cnt==nc
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
%     if sum(EE==srtrg(n))==0

%     end
    if  nop==nopM
        break
    end
end
E
nnn=sort(c(:));
nnn(1:nc)'
