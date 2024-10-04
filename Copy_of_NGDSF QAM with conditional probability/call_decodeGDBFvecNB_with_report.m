clear
save_rslt = 1;
ZERO = 1; % all zeros seq
p = 4;
q = 2^p;
pth1 = (fullfile(pwd, 'related_functions'));
addpath(pth1);
pth2 = (fullfile(pwd, 'related_variables'));
pth3 = (fullfile(pwd, 'related_variables/GF_arithm'));
pth4 = (fullfile(pwd, 'related_variables/alists'));
pth5 = (fullfile(pwd, 'related_variables/alists/matrices'));
pth6 = (fullfile(pwd, 'results/'));
words = (0:q-1);
H_matrix_mat_fl_nm = '204.102.3.6.16';
load([fullfile(pth4, H_matrix_mat_fl_nm) '.mat']);
h = full(h);
hl=double(h>0);
N = size(h,2);
M = size(h,1);
K = N-M;  %%%%%% important

Rate = K/N;

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
parf = 1200;
ebn0              = 5.6:0.2:7.6;  
snr_cnt           = length(ebn0);
max_gen           = 1e6;
max_err_cnt1      = 60; % catch max_err_cnt1 error frames at low Eb_No(<Eb_No_thrshld)
max_err_cnt2      = 60; % catch max_err_cnt1 error frames at high Eb_No
Eb_No_thrshld     = 4.8; % low-high Eb_No threshold 
mnk               = 25; % flip the mnk symbols that have the lowest flipping metric 
max_iter          = 1000;    % Max iterations for decoding
max_max_iter      = 1200; % keep try above max_iter with single VN permute until l = max_max_l
theta             = nan;   % Flipping threshold
eta               = nan;     % Perturbation noise scale parameter
eta1               = nan; % QAM conditional probability noise factor

SNR_lev = [9.6;9.8;10.2;10.4; 10.6];

conf_tbl = [0.05 0.9 1.2;...
    0.2 0.92 1.22;...
    0.3 0.95 1.27;...
    0.37 0.96 1.32;...
    0.5 0.96 1.4;...
    
    ];
%%
[norm_const, gray_labels, I, Q] = qam_constellation(q);
gray_labels = 0:q-1; % if no need for gray coding, if needed comment this line
avg_pw = norm_const'*norm_const/q; % Symbols average power 
fct1 = sqrt(avg_pw); % to use integer constellation diagram, normalise noise amplitude by this factor
%%
p1 = p; % for QAM p1 = log2(q), for BPSK p1 = 1
EbNo_linear = 10.^(ebn0/10);
symbolEnergy = 1; % in symbols are normalised to have average power of 1
noiseVariance = symbolEnergy ./ (p1 * EbNo_linear * Rate);
Snr_db = 10*log10(1./noiseVariance);
sigma0 = sqrt(noiseVariance/2) ;

sigma = fct1*sigma0; % as symbols aren't normalized, noise is multiplied by the inversse of symbol average power
%%
rng(0); % repetitive noise generation
rng(1);
alph_bin =  fliplr(dec2bin(words, p) - 48);
alph_bin_mod = (-1).^alph_bin;

NFrames(snr_cnt)=0;
if ZERO
    [G,~] = Generator_matrix_G_from_full_rank_H(h, add_mat, mul_mat, div_mat);
end
%%
clear conf_detail
conf_detail.a10 = sprintf("algoritm : %s",'Non Binary R-QAM GD Symbol Flipping');
conf_detail.a11 = sprintf("H matrix : %s",H_matrix_mat_fl_nm);
conf_detail.a12 = sprintf("N = %d, M = %d, K = %d, GF(%d)",N,M,K,q);
conf_detail.a13 = sprintf("max iter : %d",max_iter);
conf_detail.a14= sprintf("max with extra iter : %d",max_max_iter);
conf_detail.a15= sprintf("max symbols to flip per iter : %d",mnk);
conf_detail.a16= sprintf("Flipping threshold : %.3f",theta);
conf_detail.a17= sprintf("perturbation noise scaling : %.3f",eta);
conf_detail.a18= sprintf("Conditional probability sigma scaling factor : %.3f",eta1);
conf_detail.a19= sprintf("max seq generation : %d", max_gen);

current_date = datestr(now, 'yyyy_mm_dd');
current_time = datestr(now, 'HH_MM_SS');
report_fle_nme0 = strcat('Non_Binary_RQAM_GD_SymboleFlipping_', H_matrix_mat_fl_nm, '_dated_' ,current_date ,'_' ,current_time);
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
    
    msg = sprintf("SNR = %.3f dB, FER = %d/%d = %.8f,// BER = %d/%d = %.8f, aver_iter = %.3f\n",...
        Snr_db(i0), FER_cnt(i0), F_gen(i0), FERR(i0), BER_cnt(i0),B_gen(i0),...
        BERR(i0), aver_iter(i0) );
    fprintf(msg)
    
    if SNR_lev(1)<Snr_db(i0)
        i1 = find(SNR_lev<Snr_db(i0),1,'last');
    else
        i1=1;
    end
    conf_tbl1 = conf_tbl(i1,:);
    theta = conf_tbl1(1);
    eta = conf_tbl1(2);
    eta1 = conf_tbl1(3);
    csigma = sigma(i0);
    nsigma    = eta*csigma;
    csigma_1 = eta1*csigma;
    while FER_cnt(i0) < max_err_cnt
        nerrS_decd_ = 0;
        nerrB_decd_ = 0;
        nerrF_decd_ = 0;
        statstc_iter__ = zeros(1,parf);
        iters_ = 0;

        info_seq = zeros(1,K);
        code_seq = zeros(1,N);
        code_seq_comp(1,1:N) = norm_const(code_seq+1);
        parfor pp=1 :parf
            if ZERO
                [info_seq, code_seq, code_seq_comp, valid_symdrom] = ...
                    generate_and_encode_QAM(ZERO, h,G, add_mat, mul_mat,p, norm_const,...
                    gray_labels);%, [10*ones(1,102) zeros(1,K-102)]);
            else
                info_seq = zeros(1,K);
                code_seq = zeros(1,N);
                code_seq_comp = norm_const(code_seq+1);
            end

            nse = csigma*randn(size(code_seq_comp))+csigma*randn(size(code_seq_comp))*1i;
            y_cmp = code_seq_comp+nse;
            [mm, HD2_cmp] = HD_QAM(y_cmp, norm_const);
            HD2=gray_labels(mm+1);

            [failed__, seqgf__, iters__] = decodeGDSFvecNB_QAM(...
                h_arr,mul_mat, add_mat,y_cmp,hl, h, max_iter,max_max_iter,...
                str_vn_cn, str_cn_vn, dc, dv,theta,nsigma,csigma_1, mnk,...
                norm_const,gray_labels, code_seq);

            seqgf__ = gray_labels(seqgf__+1)';

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
        %
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
        msg = sprintf("SNR = %.3f dB, FER = %3d/%8d = %.3e,// BER = %4d/%11d = %.3e, aver_iter = %.3f\n",...
            Snr_db(i0), FER_cnt(i0), F_gen(i0), FERR(i0), BER_cnt(i0),B_gen(i0),...
            BERR(i0), aver_iter(i0) );
        fprintf(msg)
        msgs  = details_in_lines(Snr_db, FER_cnt,BER_cnt, SER_cnt, F_gen, K,p, aver_iter, conf_detail, report_fle_nme, save_rslt);
        if save_rslt
            save([report_fle_nme,'.mat']);
        end
    end
end