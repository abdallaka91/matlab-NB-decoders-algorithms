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
Nblocks = 1000;

for blk = 1 : Nblocks
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


    L = zeros(n+1, N,q);
    HL = zeros(n+1, N);
    L(1,:,:) = prb;


    L(2,1:4,:) = f_f_q(add_mat, sqz(L(1,:,:),1), q);
    L(3,1:2,:) = f_f_q(add_mat, sqz(L(2,1:4,:),1), q);
    L(4,1,:) = f_f_q(add_mat, sqz(L(3,1:2,:),1), q);

    HL(4,1)  = HD_N_q(sqz(L(4,1,:),1));
    HL(4,1) = 0;
    L(4,2,:) = f_g_q(HL(4,1), sqz(L(3,1:2,:),1) , q, add_mat);
    HL(4,2) = HD_N_q(sqz(L(4,2,:),1));
    HL(4,2) = 0;%
    HL(3,1:2) = f_c_q(add_mat, sqz(HL(4,1:2),1));

    L(3,3:4,:) = f_g_q(HL(3,1:2), sqz(L(2,1:4,:),1) , q, add_mat);
    L(4,3,:) = f_f_q(add_mat, sqz(L(3,3:4,:),1), q);
    HL(4,3)  = HD_N_q(sqz(L(4,3,:),1));
    HL(4,3) = 0;
    L(4,4,:) = f_g_q(HL(4,3), sqz(L(3,3:4,:),1) , q, add_mat);
    HL(4,4) = HD_N_q(sqz(L(4,4,:),1));

    HL(3,3:4) = f_c_q(add_mat, sqz(HL(4,3:4),1));
    HL(2,1:4) = f_c_q(add_mat, sqz(HL(3,1:4),1));
    %--------------------------------------------------------------

    L(2,5:8,:) = f_g_q(HL(2,1:4), sqz(L(1,:,:),1) , q, add_mat);

    L(3,5:6,:) = f_f_q(add_mat, sqz(L(2,5:8,:),1), q);
    L(4,5,:) = f_f_q(add_mat, sqz(L(3,5:6,:),1), q);
    HL(4,5)  = HD_N_q(sqz(L(4,5,:),1));
    HL(4,5) = 0;
    L(4,6,:) = f_g_q(HL(4,5), sqz(L(3,5:6,:),1) , q, add_mat);
    HL(4,6) = HD_N_q(sqz(L(4,6,:),1));

    HL(3,5:6) = f_c_q(add_mat, sqz(HL(4,5:6),1));

    L(3,7:8,:) = f_g_q(HL(3,5:6), sqz(L(2,5:8,:),1) , q, add_mat);
    L(4,7,:) = f_f_q(add_mat, sqz(L(3,7:8,:),1), q);
    HL(4,7)  = HD_N_q(sqz(L(4,7,:),1));

    L(4,8,:) = f_g_q(HL(4,7), sqz(L(3,7:8,:),1) , q, add_mat);
    HL(4,8) = HD_N_q(sqz(L(4,8,:),1));
    %--------------------------------------------------------------
    dcd = HL(4,:);
    msg_cap = dcd(Q1(N-K+1:end));
    Nerrs = sum(msg~=msg_cap);

    if Nerrs>0
        Nbiterrs = Nbiterrs +  Nerrs;
        Nblkerrs = Nblkerrs + 1;
    end
end

BER_sim = Nbiterrs/K/Nblocks
FER_sim = Nblkerrs/Nblocks
dcd1(Q1) = dcd;
msg_rec= dcd1(end-K+1:end);
msg - msg_rec


