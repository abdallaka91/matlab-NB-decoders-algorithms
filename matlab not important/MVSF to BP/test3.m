clear;
q = 64;
pth1 = (fullfile(pwd, 'related_functions'));
addpath(pth1);
pth2 = (fullfile(pwd, 'related_variables'));
pth3 = (fullfile(pwd, 'related_variables/GF_arithm'));
pth4 = (fullfile(pwd, 'related_variables/alists'));
pth5 = (fullfile(pwd, 'related_variables/alists/matrices'));
pth6 = (fullfile(pwd, 'results/'));
H_matrix_mat_fl_nm = 'BeiDou_44_bb_GF64';

load([fullfile(pth4, H_matrix_mat_fl_nm) '.mat']);

fl_nm = ['arith_' num2str(q) '.mat'];
if  exist(fullfile(pth3, fl_nm), 'file') == 2
    load(fullfile(pth3, fl_nm));
else
    add_mat = GF_arithm_matrix(q, 'add');
    mul_mat = GF_arithm_matrix(q, 'mul');
    div_mat = GF_arithm_matrix(q, 'div');
    save(fullfile(pth3, ['arith_' num2str(q) '.mat']), 'add_mat' ,'mul_mat','div_mat')
end

words = (0:q-1);
p=log2(q);
ZERO=1; 
h = full(h);
N = size(h,2);
M = size(h,1);
 K = N-M;
 nm = 16;
nc2v=nm;
cmp_c2v =150;
dc1 = [0 4];
Dev_pos_cnt = length(dc1)-1;
di = cell(length(dc1),1);
di{1} = [0 0];
di{2} = [0 15 2 2 1];
% di{3} = [0 2 2];
v_weights = [1.5 0.5]*256;
LLRfact = 1024;
unreliable_sat=-inf;
parforN = 50;
max_err_cnt1 = 60;
max_err_cnt2 = 30; 
Eb_No_thrshld = 3.8;
max_gen = 1e6;
max_iter = 50;
max_attempt = 1;
ebn0 = 1.6; 

dev_lsts = cell(M,1);
dev_pos = cell(M,1);
str_cn_vn = cell(M,1);
dc = zeros(M,1);
for i = 1 : M
    str_cn_vn{i, 1} = find(h(i,:));
    dc(i) = length(str_cn_vn{i});
    dc11 = dc1;
    dc11(1) = dc(i)-sum(dc1);
    [lst_deviation_lst, lst_dev_pos, dev_lsts_i, dev_pos_i]=list_dev_reliabl(dc11, di);
    dev_lsts{i} = dev_lsts_i;
    dev_pos{i} = dev_pos_i;
end


p1 = 1;
Rate = p1*K/N; 
ebn0_n = 10.^(ebn0/10);
N0 = 1./(Rate*ebn0_n);
sigma = sqrt(N0/2);
snr = -10*log10(2*sigma.^2);

alph_bin =  fliplr(dec2bin(words, p) - 48);
alph_bin_mod = (-1).^alph_bin;


str_vn_cn = cell(N,1);
dv = zeros(N,1);
for j = 1 : N
    str_vn_cn{j, 1} = (find(h(:,j)))';
    dv(j) = length(str_vn_cn{j});
end
snr_cnt = length(sigma);
FERstat = zeros(snr_cnt,1);
SERstat = zeros(snr_cnt,1);
BERstat = zeros(snr_cnt,1);
gen_seq_cnt = 0*ones(snr_cnt,1);
FER = zeros(snr_cnt,1);
SER = zeros(snr_cnt,1);
BER = zeros(snr_cnt,1);
iter_cnt = zeros(snr_cnt,1);
aver_iter = zeros(snr_cnt,1);
max_err_cnt = max_err_cnt1;

needed_iters = nan(max_gen,snr_cnt);
LLR_20 =zeros(N,q);

[G,~] = Generator_matrix_G_from_full_rank_H(h, add_mat, mul_mat, div_mat);
rng(0)
for i0 = 1 : snr_cnt
    if ebn0(i0)>=Eb_No_thrshld
        max_err_cnt = max_err_cnt2;
    end
    iter_cnt_ = 0;
    FER_ = 0;
    SER_ = 0;
    BER_ = 0;
    gen_seq_cnt_ = 0;
    needed_iters_ = nan(parforN,1);
    msg = sprintf("EbNo = %.3f dB, FER = %d/%d = %.8f,// BER = %d/%d = %.8f, aver_iter = %.3f\n",...
        ebn0(i0), FER(i0), gen_seq_cnt(i0), FER(i0)/gen_seq_cnt(i0), BER(i0), gen_seq_cnt(i0)*K*p,...
        FER(i0)/(gen_seq_cnt(i0)*K*p), 0);
    fprintf(msg)
    sigm =sigma(i0);
    KK=0;
    while FER(i0) < max_err_cnt && gen_seq_cnt(i0)<max_gen
        parfor j = 1 : parforN
            if ZERO
                [info_seq, code_seq, valid_symdrom, y_bin] = generate_and_encode(ZERO, h,G, add_mat, mul_mat, p);
            else
                info_seq = zeros(1,K);
                code_seq = zeros(1,N);
                y_bin = ones(N,p);
            end
            gen_seq_cnt_ = gen_seq_cnt_+1;
            nse = sigm*randn(size(y_bin));
            y_bin_nse = y_bin + nse;
% load y_bin_nse.mat
            LLR_2 = -LLR_simple3(y_bin_nse,LLRfact , unreliable_sat, q,N, alph_bin, LLR_20);
            [~,HD1] = min(LLR_2,[], 2);
            HD1 = HD1'-1;
            nes = sum(HD1~=code_seq);
            [iters, dec_seq, success_dec] = ...
                presorted_MVSF_with_rinfrc(LLR_2, max_iter, mul_mat, add_mat, div_mat,...
                h,str_cn_vn, dc, dev_lsts,dv, str_vn_cn, nm, v_weights, max_attempt, cmp_c2v, nc2v);
            needed_iters_(j) = iters;

            iter_cnt_ = iter_cnt_ + iters;
            rec_info_seq = dec_seq(1:K);
            nd = sum(rec_info_seq~=info_seq);
            if ~success_dec && nd ~=0
                FER_ = FER_ +1;
                [SER0, BER0] = SER_BER(info_seq,rec_info_seq,p, 0, 0);
                SER_ = SER0 + SER_;
                BER_ = BER0 + BER_;
            end

        end
        needed_iters(KK+1:KK+parforN,i0) = needed_iters_;
        KK=KK+parforN;
        SER(i0) = SER_;
        BER(i0) = BER_;
        gen_seq_cnt(i0) = gen_seq_cnt_;
        iter_cnt(i0) = iter_cnt_;
        aver_iter(i0) = iter_cnt(i0)/gen_seq_cnt(i0);

        FER(i0) = FER_;
        FERstat(i0)=FER(i0)/gen_seq_cnt(i0);
        SERstat(i0)=SER(i0)/(gen_seq_cnt(i0)*K);
        BERstat(i0)=BER(i0)/(gen_seq_cnt(i0)*K*p);

        fprintf(repmat('\b',1,length(char(msg))));
        msg = sprintf("EbNo = %.3f dB, FER = %d/%d = %.8f,// BER = %d/%d = %.8f, aver_iter = %.3f\n",...
            ebn0(i0), FER(i0), gen_seq_cnt(i0), FER(i0)/gen_seq_cnt(i0), BER(i0), gen_seq_cnt(i0)*K*p,...
            FER(i0)/(gen_seq_cnt(i0)*K*p), aver_iter(i0) );
        fprintf(msg)
    end
    
end