clear;%/home/abdallah/Downloads/NB_LDPC_Decoders/call_MV_SF_parforloop.m
comput_SER_BER = false;
ZERO=1; % if  0 then simulate all zeros sequence
plt = 0; % continiously plot FER performance if 1
nm = 4;% V2C m2ssage size
dc1 = [0 1 2];
save_rslt = 0;
Dev_pos_cnt = length(dc1)-1;
di = cell(length(dc1),1);
di{1} = [0 0];
di{2} = [0 1];
di{3} = [0 2 1];
v_weights = [0.8 0.3]*256; %good forDV=3
parforN = 100;

% di{1} = [0 0];
% di{2} = [0 1];
% di{3} = [0 2 1];
% v_weights = [140 160];%good forDV=2

max_err_cnt1 = 60; % at low Eb_No(<Eb_No_thrshld)
max_err_cnt2 = 30; %at high Eb_No

Eb_No_thrshld = 3.8;
max_gen = 1e6;
max_iter = 25;
ebn0 = 3:0.2:4.8; %dB
p = 4;
q = 2^p;
PERM_rng = 80;
PERM = 40;
max_trial = 10;

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
% h=H;
h = full(h);
N = size(h,2);
M = size(h,1);
K = N-M;
fl_nm = ['arith_' num2str(q) '.mat'];
if  exist(fullfile(pth3, fl_nm), 'file') == 2
    load(fullfile(pth3, fl_nm));
else
    add_mat = GF_arithm_matrix(q, 'add');
    mul_mat = GF_arithm_matrix(q, 'mul');
    div_mat = GF_arithm_matrix(q, 'div');
    save(fullfile(pth3, ['arith_' num2str(q) '.mat']), 'add_mat' ,'mul_mat','div_mat')
end

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
disp(length(dev_lsts_i))

%%
clear conf_detail
conf_detail.a11fl_nme = sprintf("H matrix : %s",H_matrix_mat_fl_nm);
conf_detail.a12Code = sprintf("N = %d, M = %d, K = %d, GF(%d)",N,M,K,q);
conf_detail.a13votesW = sprintf('V0 = %.3f, V1 = %.3f', v_weights(1) , v_weights(2));
conf_detail.a14iter = sprintf("max iter : %d",max_iter);
conf_detail.a15max_trail = sprintf("max nb of trials : %d",max_trial);
conf_detail.a16max_seq = sprintf("max seq generation : %d", max_gen);
conf_detail.a17fl_nme = sprintf("max error frame detection: %d",max_err_cnt1);
conf_detail.a18dev_cnt = sprintf("nb of deviation paths : %d", length(dev_lsts_i));
conf_detail.a19reg_widths = sprintf("regions width (high to low reliable) : {") + sprintf(" %d ", dc11) + sprintf("}");
formatted_str = '';
for i = 1:length(di)
    formatted_str = [formatted_str  sprintf(['[ ' repmat('%d ',1,  length(di{i})-1) ']'], di{i}(2:end))];
    if i < length(di)
        formatted_str = [formatted_str ' || '];
    end
end
conf_detail.a21di1_di2 = sprintf("%s\n", formatted_str);

current_date = datestr(now, 'yyyy_mm_dd');
current_time = datestr(now, 'HH_MM_SS');
report_fle_nme0 = strcat(extractAfter(conf_detail.a11fl_nme, "H matrix : "), '_' ,current_date ,'_' ,current_time);
report_fle_nme = pth6+report_fle_nme0;
%%

p1 = 1;
Rate = p1*K/N; %p1 is nb of bits per channel use with the modulation, for example for bpsk it is 1

ebn0_n = 10.^(ebn0/10);
N0 = 1./(Rate*ebn0_n);
sigma = sqrt(N0/2);
snr = -10*log10(2*sigma.^2);

%%
alph_bin =  fliplr(dec2bin(words, p) - 48);
alph_bin_mod = (-1).^alph_bin;
[G,~] = Generator_matrix_G_from_full_rank_H(h, add_mat, mul_mat, div_mat);
% info_seq = [12,5,15,5,6,14,15,11,1,14,5,12,14,7,15,9,8,9,6,12,1,3,15,9,13,11,8,15,1,8,1,3,2,4,8,13,4,11,10,12,11,11,2,11,3,2,9,2,9,5,6,9,5,8,3,7,13,1,9,13,3,8,6,10,3,12,0,8,1,14,4,15,2,0,12,10,12,9,11,11,15,10,12,10,15,12,14,10,10,0,8,10,4,8,12,14,8,8,11,1,6,12];
% code_seq = gf_mat_mul(info_seq,G, add_mat, mul_mat);
% valid_symdrom = gf_mat_mul(code_seq,h', add_mat, mul_mat);
% y_bin0 = fliplr(dec2bin(code_seq, p) - 48);
% y_bin = (-1).^y_bin0;

%%

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

for i0 = 1 : snr_cnt
    if ebn0(i0)>=Eb_No_thrshld
        max_err_cnt = max_err_cnt2;
    end
    iter_cnt_ = 0;
    FER_ = 0;
    FER__ = 0;
    SER_ = 0;
    BER_ = 0;
    gen_seq_cnt_ = 0;
    msg = sprintf("EbNo = %.3f dB, FER = %d/%d = %.8f,// BER = %d/%d = %.8f, aver_iter = %.3f\n",...
        ebn0(i0), FER(i0), gen_seq_cnt(i0), FER(i0)/gen_seq_cnt(i0), BER(i0), gen_seq_cnt(i0)*K*p,...
        FER(i0)/(gen_seq_cnt(i0)*K*p), 0);
    fprintf(msg)
    sigm =sigma(i0);
    while FER(i0) < max_err_cnt && gen_seq_cnt(i0)<max_gen
        parfor j = 1 : parforN
            [info_seq, code_seq, valid_symdrom, y_bin] = generate_and_encode(ZERO, h,G, add_mat, mul_mat, p);
            gen_seq_cnt_ = gen_seq_cnt_+1;
            nse = sigm*randn(size(y_bin));
            y_bin_nse = y_bin + nse;
            LLRfact = 1;
            unreliable_sat=-inf;
            LLR_2 = 1024*LLR_simple3(y_bin_nse, p,LLRfact , unreliable_sat);
            LLR_2 = round(LLR_2);
            [iters,Trial, dec_seq, success_dec] = presorted_MVSF_3(LLR_2, max_iter, mul_mat, add_mat, div_mat, h,str_cn_vn, dc, ...
                dev_lsts,dev_pos,nm,v_weights, Dev_pos_cnt, PERM_rng, PERM, max_trial);
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

        msgs  = details_in_lines(ebn0, FER,BER, SER, gen_seq_cnt, K,p, aver_iter, conf_detail, report_fle_nme, save_rslt);
        if save_rslt
            save(report_fle_nme+'.mat');
        end

    end
end
disp(msgs{end})