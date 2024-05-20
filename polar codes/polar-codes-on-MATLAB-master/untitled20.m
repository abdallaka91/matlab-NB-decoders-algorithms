clear
close all



snr=6;

clear
load Q_reliability.mat
N = 16;
K = 8;
EbNodB = 4;
Rate = K/N;
EbNo=10^(EbNodB/10);
sigma=sqrt(1/(2*Rate*EbNo));
n = log2(N);
Q1 = Q(Q<=N);
M=N-K;
AC= sort(Q1(1:M));
A=setdiff(1:N,AC);
Nbiterrs = 0;
Nblkerrs = 0;
Nblocks = 10000;
u_A=randi([0 1],1,K);  % information vetor
u_AC=zeros(1,K);  % frozen vector

[u,x] = polar_code_encoder(n,A,u_A,AC,u_AC);
y = polar_code_channel(N,x,inf);
y=y+sigma*randn(1,N);
u_e = polar_code_SC_decoder(n,N,y,AC);

if(u==u_e)
    disp('correct');
else
    disp('error');
end