clear
load Q_reliability.mat
N = 4;
K = 2;
EbNodB = -1;
Rate = K/N;
EbNo=10^(EbNodB/10);
sigma=sqrt(1/(2*Rate*EbNo));
n = log2(N);
Q1 = Q(Q<=N);
F= Q1(1:N-K);

Nbiterrs = 0;
Nblkerrs = 0;
Nblocks = 1000;

for blk = 1 : Nblocks


    msg =randi ([0 1],1,K);

    u = zeros(1, N);
    u(Q1(N-K+1:end)) = msg;
    m = 1;
    for d=n-1:-1:0
        for i=1:2*m:N
            a = u(i:i+m-1);
            b= u(i+m:i+2*m-1);
            u(i:i+2*m-1)=[mod(a+b, 2) b]; %combining

        end
        m=m*2;
    end
    cword = u;
    s = 1-2*cword;
    r=s+sigma*randn(1,N);

    g = @(a,b,c) b+(1-2*c).*a;
    f = @(a,b) (1-2*(a<0)).*(1-2*(b<0)).*min(abs(a), abs(b));
    f1 = @(a,b) [xor(a, b) b];

    L1l = [f(r(1),r(3)), f(r(2), r(4))];
    L2ll=f(L1l(1), L1l(2));
    ur(1) = double(L2ll<0);
    ur(1) = 0; %frozen

    L2lr=g(L1l(1), L1l(2), ur(1));
    ur(2) = double(L2lr<0);
    ur(2) = 0; %frozen

    v1l = [xor(ur(1), ur(2)) ur(2)];

    L1r = [g(r(1),r(3), v1l(1)), g(r(2), r(4), v1l(2))];
    L2ll=f(L1r(1), L1r(2));
    ur(3) = double(L2ll<0);

    L2lr=g(L1r(1), L1r(2),ur(3));
    ur(4) = double(L2lr<0);

    v1r = [xor(ur(3), ur(4)) ur(4)];

    ucap = [xor(v1l, v1r), v1r];

    msg_cap = ur(Q1(N-K+1:end));
    Nerrs = sum(msg~=msg_cap);

    if Nerrs>0
        Nbiterrs = Nbiterrs +  Nerrs;
        Nblkerrs = Nblkerrs + 1;
    end
end

BER_sim = Nbiterrs/K/Nblocks
FER_sim = Nblkerrs/Nblocks
