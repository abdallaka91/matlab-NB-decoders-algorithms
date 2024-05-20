clear
pth1 = (fullfile(pwd, 'related_functions'));
addpath(pth1);
pth2 = (fullfile(pwd, 'related_variables'));
pth3 = (fullfile(pwd, 'related_variables/GF_arithm'));

q=4;
p = log2(q);
words = (0:q-1);

fl_nm = ['arith_' num2str(q) '.mat'];
if  exist(fullfile(pth3, fl_nm), 'file') == 2
    load(fullfile(pth3, fl_nm));
else
    add_mat = GF_arithm_matrix(q, 'add');
    mul_mat = GF_arithm_matrix(q, 'mul');
    div_mat = GF_arithm_matrix(q, 'div');
    save(fullfile(pth3, ['arith_' num2str(q) '.mat']), 'add_mat' ,'mul_mat','div_mat')
end

load Q_reliability.mat
N = 8;
K = 4;
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
Nblocks = 10000;

msg=[1 ];
[cword, msg] = gen_rand_nb_code(N,K, Q1, q, add_mat, mul_mat,msg);


y_bin0 = fliplr(dec2bin(cword, p) - 48);
y_bin = (-1).^y_bin0;
nse = sigma*randn(size(y_bin));
y_bin_nse = y_bin + nse;

prb = prb_q(y_bin_nse, sigma);
% LLR_2 = LLR_simple3(y_bin_nse, p,1 , -1000);
[mm,HD_seq] = max(prb,[],2);
HD_seq = HD_seq'-1;
ndff = sum(HD_seq~=cword);


decw = polar_nb_dec(add_mat, prb,frozen_idx,M,N,q,n );

msg_cap = decw(Q1(N-K+1:end));
Nerrs = sum(msg~=msg_cap);

if Nerrs>0
    Nbiterrs = Nbiterrs +  Nerrs;
    Nblkerrs = Nblkerrs + 1;
end