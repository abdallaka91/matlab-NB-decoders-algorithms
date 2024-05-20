clear
load Q_reliability.mat
N = 64;
K = 32;
EbNodB = 3;
Rate = K/N;
EbNo=10^(EbNodB/10);
sigma=sqrt(1/(2*Rate*EbNo));
n = log2(N);
Q1 = Q(Q<=N);
M=N-K;
frozen_idx= sort(Q1(1:M));
Nbiterrs = 0;
Nblkerrs = 0;
Nblocks = 1000;


stat_v = cell(n+1,1);
L=stat_v;
V = L;
for i = 1 : n+1
    i0 = i-1;
    L{i,1} = nan(2^i0,  2^(n-i0));
    V{i,1} = nan(2^i0,  2^(n-i0));
    stat_v{i,1} = false(1,  2^(i0));
end


for blk = 1 : Nblocks
    [cword, msg] = gen_rand_code(N,K, Q1);
    s = 1-2*cword;
    rn=s+sigma*randn(1,N);
    % L{1} = rn;
    % stat_L{1} = true;
    %
    % i = 1;
    % j = 1;
    % froz_j1 = 1;
    % decw = nan(1,16);
    % while i>0
    %     if stat_v{i}(j)
    %         if mod(j,2)==0
    %             i=i-1;
    %             j=j/2;
    %             temp0 = V{i+1}(2*j-1,:);
    %             temp1 = V{i+1}(2*j,:);
    %             temp = function_comb(temp0, temp1);
    %             V{i}(j,:) = temp;
    %             stat_v{i}(j) = true;
    %
    %         else
    %             i = i-1;
    %             j=(j+1)/2;
    %         end
    %     elseif stat_v{i+1}(2*j-1)
    %         i = i+1;
    %         j=2*j;
    %         temp0 = L{i-1}(j/2,:);
    %         temp1 = V{i}(j-1,:);
    %         temp = function_g(temp0, temp1);
    %         L{i}(j,:) = temp;
    %         if i==n+1
    %             stat_v{i}(j)=true;
    %             Li = L{i}(j,:);
    %             V{i}(j,1) = double(Li<0);
    %             if froz_j1<=M
    %                 if j==frozen_idx(froz_j1)
    %                     V{i}(j,1) = 0;
    %                     froz_j1 = froz_j1+1;
    %                 end
    %             end
    %             decw(1:j) = V{i}(1:j);
    %         end
    %     else
    %         i=i+1;
    %         j=2*j-1;
    %         temp0 = L{i-1}((j+1)/2,:);
    %         temp = function_f(temp0);
    %         L{i}(j,:) = temp;
    %         if i==n+1
    %             stat_v{i}(j)=true;
    %             Li = L{i}(j,:);
    %             V{i}(j,1) = double(Li<0);
    %             if froz_j1<=M
    %                 if j==frozen_idx(froz_j1)
    %                     V{i}(j,1) = 0;
    %                     froz_j1 = froz_j1+1;
    %                 end
    %             end
    %             decw(1:j) = V{i}(1:j);
    %             i=i-1;
    %             j=(j+1)/2;
    %         end
    %     end
    % end
    decw = polar_dec( rn, L, stat_v, V,frozen_idx,M,n );

    msg_cap = decw(Q1(N-K+1:end));
    Nerrs = sum(msg~=msg_cap);

    if Nerrs>0
        Nbiterrs = Nbiterrs +  Nerrs;
        Nblkerrs = Nblkerrs + 1;
    end
end

BER_sim = Nbiterrs/K/Nblocks;
FER_sim = Nblkerrs/Nblocks;



