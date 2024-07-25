clear
save_rslt = 1;
ZERO = 1; % all zeros seq
p = 6;
q = 2^p;
pth1 = (fullfile(pwd, 'related_functions'));
addpath(pth1);
pth2 = (fullfile(pwd, 'related_variables'));
pth3 = (fullfile(pwd, 'related_variables/GF_arithm'));
pth4 = (fullfile(pwd, 'related_variables/alists'));
pth5 = (fullfile(pwd, 'related_variables/alists/matrices'));
pth6 = (fullfile(pwd, 'results/'));
words = (0:q-1);
H_matrix_mat_fl_nm = '384.320.4.20.64';
load([fullfile(pth4, H_matrix_mat_fl_nm) '.mat']);
h = full(h);
hl=double(h>0);
N = size(h,2);
M = size(h,1);
K = N-M;  %%%%%% important


Nb = N*p;
Kb = K*p;
fl_nm = ['arith_' num2str(q) '.mat'];
if  exist(fullfile(pth3, fl_nm), 'file') == 2
    load(fullfile(pth3, fl_nm));
else
    add_mat = GF_arithm_matrix(q, 'add');
    mul_mat = GF_arithm_matrix(q, 'mul');
    div_mat = GF_arithm_matrix(q, 'div');
    save(fullfile(pth3, ['arith_' num2str(q) '.mat']), 'add_mat' ,'mul_mat','div_mat')
end

dv = zeros(1,N);
str_vn_cn = cell(N,1);
for j = 1 : N
    str_vn_cn{j, 1} = find(h(:,j));
    dv(j) = length(str_vn_cn{j});
end

str_cn_vn = cell(M,1);
dc = zeros(M,1);
h_sz = 0;
for i = 1 : M
    str_cn_vn{i, 1} = find(h(i,:));
    dc(i) = length(str_cn_vn{i});
    h_sz = h_sz + dc(i) ;
end

h_arr =zeros(1,h_sz*3);
k=0;
for i = 1 : M
    for j = 1 : dc(i)
        h_arr(k*3+1:(k+1)*3) =[i str_cn_vn{i}(j) h(i,str_cn_vn{i}(j))];
        k=k+1;
    end
end
parf = 40;
ebn0              = 4:0.2:5.6;  % SNR values to simulate
snr_cnt           = length(ebn0);
max_gen           = 1e6;
max_err_cnt1      = 50; % at low Eb_No(<Eb_No_thrshld)
max_err_cnt2      = 50; %at high Eb_No
Eb_No_thrshld     = 4.8;
mnk               = 15;
max_iter          = 1000;    % Max iterations for decoding
max_max_iter      = 2000; % keep try above T with single VN permute until l = max_max_l
syndrome_weight   = 0.95;   % Syndrome weight parameter
theta             = -2.7;   % Flipping threshold
eta               = 1.1;     % Perturbation noise scale parameter
nones1            = 4; %max hamming dist
nones2            = 1; %max hamming dist
sngl_VN_th        = 4.8; %dB above whichSNR check if a single symbol flipping could reduce nb of faild by its dv


GDBFerrors       = zeros(snr_cnt, 1);
GDBFframeerrors  = zeros(snr_cnt, 1);
GDBFframeundetected  = zeros(snr_cnt, 1);
NFrames          = zeros(snr_cnt, 1);

p1 = 1;
Rate = p1*K/N; %p1 is nb of bits per channel use with the modulation, for example for bpsk it is 1
ebn0_n = 10.^(ebn0/10);
N0 = 1./(Rate*ebn0_n);
sigma = sqrt(N0/2);

%%
rng(1); % repetitive noise generation
alph_bin =  fliplr(dec2bin(words, p) - 48);
alph_bin_mod = (-1).^alph_bin;

NFrames(snr_cnt)=0;
if ZERO
    [G,~] = Generator_matrix_G_from_full_rank_H(h, add_mat, mul_mat, div_mat);
end
dev1 =     deviation_hamm_lst(p,nones1);
dev_hamm1 = dev1==1;
dev2 =     deviation_hamm_lst(p,nones2);
dev_hamm2 = dev2==1;

%%
clear conf_detail
conf_detail.a10 = sprintf("algoritm : %s",'Non Binary BPSK GDBF');
conf_detail.a11 = sprintf("H matrix : %s",H_matrix_mat_fl_nm);
conf_detail.a12 = sprintf("N = %d, M = %d, K = %d, GF(%d)",N,M,K,q);
conf_detail.a13 = sprintf("max iter : %d",max_iter);
conf_detail.a14= sprintf("max with extra iter : %d",max_max_iter);
conf_detail.a15= sprintf("max bits to flip per iter : %d",mnk);
conf_detail.a16= sprintf("Binary syndrom weight : %.3f",syndrome_weight);
conf_detail.a17= sprintf("Flipping threshold : %.3f",theta);
conf_detail.a18= sprintf("perturbation noise scaling : %.3f",eta);
conf_detail.a19= sprintf("max seq generation : %e", max_gen);
conf_detail.a20= sprintf("single VN symbol flipping max Hamming distance 1: %d", nones1);
conf_detail.a21= sprintf("single VN symbol flipping max Hamming distance 2: %d", nones2);
conf_detail.a22= sprintf("single VN symbol flipping applying threshold: %d", sngl_VN_th);



current_date = datestr(now, 'yyyy_mm_dd');
current_time = datestr(now, 'HH_MM_SS');
report_fle_nme0 = strcat('Non_Binary_BPSK_GDBF_', H_matrix_mat_fl_nm, '_dated_' ,current_date ,'_' ,current_time);
report_fle_nme = fullfile(pth6, report_fle_nme0);
%%


SER_cnt=zeros(snr_cnt,1);
BER_cnt=zeros(snr_cnt,1);
FER_cnt=zeros(snr_cnt,1);

S_gen=zeros(snr_cnt,1);
B_gen=zeros(snr_cnt,1);
F_gen=zeros(snr_cnt,1);

SERR=zeros(snr_cnt,1);
BERR=zeros(snr_cnt,1);
FERR=zeros(snr_cnt,1);

aver_iter = zeros(snr_cnt,1);

statstc_iter = zeros(snr_cnt,max_max_iter);

iter_cnt = zeros(snr_cnt,1);
pw = (2.^(0:p-1))';

for i0=1 : snr_cnt
    max_err_cnt = max_err_cnt1;
    if ebn0(i0)>=Eb_No_thrshld
        max_err_cnt = max_err_cnt2;
    end
    csigma = sigma(i0);
    nsigma    = eta*csigma;
    msg = sprintf("EbNo = %.3f dB, FER = %d/%d = %.8f,// BER = %d/%d = %.8f, aver_iter = %.3f\n",...
        ebn0(i0), FER_cnt(i0), F_gen(i0), FERR(i0), BER_cnt(i0),B_gen(i0),...
        BERR(i0), aver_iter(i0) );

    fprintf(msg)
    keepGoing=1;
    if csigma >= sngl_VN_th
        sngl_VN = 1;
    else
        sngl_VN = 0;
    end
    while FER_cnt(i0) < max_err_cnt

        nerrS_decd_ = 0;
        nerrB_decd_ = 0;
        nerrF_decd_ = 0;
        statstc_iter__ = zeros(1,parf);
        iters_ = 0;

        parfor pp=1 :parf
            if ZERO
                [info_seq, code_seq, valid_symdrom, y_bin_mod2D] = generate_and_encode(ZERO, h,G, add_mat, mul_mat, p);
            else
                info_seq = zeros(1,K);
                code_seq = zeros(1,N);
                y_bin_mod2D = ones(N,p);
            end
            y_bin_mod1D = reshape(y_bin_mod2D', 1, []);
            y_bin2D = (1-y_bin_mod2D)/2;
            y_bin1D = reshape(y_bin2D', 1, []);
            noisevec = csigma*randn(1,Nb);
            y        = y_bin_mod1D+ noisevec;
            % load yy.mat

            [~, failed__, ~, ~, E__, ~, seqgf__, iters__] = decodeGDBFvecNB(p,pw, h_arr,str_cn_vn,dc,dv,...
                mul_mat, add_mat, y,hl, h, N, max_iter, syndrome_weight, theta,mnk,nsigma,dev_hamm1, dev_hamm2, sngl_VN_th, max_max_iter);

            statstc_iter__(1,pp) = iters__;
            rec_info_seq = seqgf__(1:K);


            if failed__~=0
                [nerrS_decd__, nerrB_decd__] = SER_BER(info_seq,rec_info_seq,p, 0, 0);
            else
                nerrS_decd__=0;
                nerrB_decd__=0;
            end
            nerrB_decd_ = nerrB_decd_ + nerrB_decd__;
            nerrS_decd_ = nerrS_decd_ + nerrS_decd__;
            iters_ = iters_ + iters__;
            if nerrB_decd__>0
                nerrF_decd_ = nerrF_decd_+1;
            end
        end
        for pp =1 : parf
            statstc_iter(i0, statstc_iter__(pp)) = statstc_iter(i0, statstc_iter__(pp)) + 1;
        end

        SER_cnt(i0) = SER_cnt(i0)+nerrS_decd_;
        BER_cnt(i0) = BER_cnt(i0)+nerrB_decd_;
        FER_cnt(i0) = FER_cnt(i0)+nerrF_decd_;

        S_gen(i0) = S_gen(i0)+K*parf;
        B_gen(i0) = B_gen(i0)+Kb*parf;
        F_gen(i0) = F_gen(i0)+parf;

        SERR(i0) = SER_cnt(i0)/S_gen(i0) ;
        BERR(i0) = BER_cnt(i0)/B_gen(i0) ;
        FERR(i0) = FER_cnt(i0)/F_gen(i0) ;

        iter_cnt(i0) = iter_cnt(i0) + iters_;
        aver_iter(i0) = iter_cnt(i0)/F_gen(i0);

        fprintf(repmat('\b',1,length(char(msg))));
        msg = sprintf("EbNo = %.3f dB, FER = %3d/%8d = %.3e,// BER = %4d/%11d = %.3e, aver_iter = %.3f\n",...
            ebn0(i0), FER_cnt(i0), F_gen(i0), FERR(i0), BER_cnt(i0),B_gen(i0),...
            BERR(i0), aver_iter(i0) );
        fprintf(msg)
        msgs  = details_in_lines(ebn0, FER_cnt,BER_cnt, SER_cnt, F_gen, K,p, aver_iter, conf_detail, report_fle_nme, save_rslt);
        if save_rslt
            save([report_fle_nme,'.mat']);
        end
    end
end