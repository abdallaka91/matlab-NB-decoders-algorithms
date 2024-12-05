clear;%/home/abdallah/Downloads/NB_LDPC_Decoders/call_MV_SF_parforloop.m
rng(1)
comput_SER_BER = false;
ZERO=0; % if  0 then simulate all zeros sequence
plt = 0; % continiously plot FER performance if 1
nm = 4;% V2C m2ssage size
p = 4;
q = 2^p;
dc1 = [0 1 2];
save_rslt = 1;
rng(1); % noise reprodudev_lstscity
Dev_pos_cnt = length(dc1)-1;
di = cell(length(dc1),1);
di{1} = [0 0];
di{2} = [0 1 ];
di{3} = [0 2 1];
v_weights = [204.8 76.8];
LLRfact = 1024;
unreliable_sat=-inf;
parforN =240;
max_err_cnt1 = 50; % at low Eb_No(<Eb_No_thrshld)
max_err_cnt2 = 50; %at high Eb_No
Eb_No_thrshld = 3.2;
max_gen = 2e5;
max_iter = 16;
max_attempt = 4;
max_max_iter=max_iter*max_attempt;
ebn0 = 1.8:0.2:3.8; %dB



projectPath = pwd; 
mainPath = fileparts(projectPath);
related_variables_pth = fullfile(mainPath, 'related_variables');
pth1 = (fullfile(pwd, 'related_functions'));
addpath(pth1);
pth3 = (fullfile(related_variables_pth, 'GF_arithm'));
pth4 = (fullfile(related_variables_pth, 'alists'));
pth5 = (fullfile(related_variables_pth, 'alists/matrices'));
pth7 = (fullfile(related_variables_pth, 'generator_matrices'));
pth6 = (fullfile(pwd, 'results/'));


words = (0:q-1);
H_matrix_mat_fl_nm = '204.102.3.6.16';
load([fullfile(pth4, H_matrix_mat_fl_nm) '.mat']);
h = full(h);
N = size(h,2);
M = size(h,1);
K = N-M;
% K=723;
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
%%
clear conf_detail
conf_detail.a11fl_nme = sprintf("H matrix : %s",H_matrix_mat_fl_nm);
conf_detail.a12Code = sprintf("N = %d, M = %d, K = %d, GF(%d)",N,M,K,q);
conf_detail.a13votesW = sprintf('V0 = %.3f, V1 = %.3f', v_weights(1) , v_weights(2));
conf_detail.a14llrfact = sprintf('LLR factor = %.3f', LLRfact);
conf_detail.a15iter = sprintf("max iter : %d",max_iter);
conf_detail.a16max_attempt = sprintf("max nb of max_attempts : %d",max_attempt);
conf_detail.a17max_seq = sprintf("max seq generation : %d", max_gen);
conf_detail.a18fl_nme = sprintf("max error frame detection: %d",max_err_cnt1);
conf_detail.a19dev_cnt = sprintf("nb of deviation paths : %d", length(dev_lsts_i));
conf_detail.a20reg_widths = sprintf("regions width (high to low reliable) : {") + sprintf(" %d ", dc11) + sprintf("}");
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
% load([related_variables_pth '/837_723_parameters_proof.mat']);
% G=double(G);
% G(:,swap)=G;
[G,~] = Generator_matrix_G_from_full_rank_H(h, add_mat, mul_mat, div_mat);
info_seq = randi([0 q-1],1,K);
code_seq = gf_mat_mul(info_seq,G, add_mat, mul_mat);
valid_symdrom = gf_mat_mul(code_seq,h', add_mat, mul_mat);
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

h0 = 0*h;
needed_iters = nan(max_gen,snr_cnt);
LLR_20 =zeros(N,q);

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
            LLR_2 = LLR_simple3(y_bin_nse,LLRfact , unreliable_sat, q,N, alph_bin, LLR_20);
            [~,HD1] = max(LLR_2,[], 2);
%             HD1 = HD1'-1;
%             nes = sum(HD1~=code_seq);

            [iters, dec_seq, success_dec] = ...
                presorted_MVSF_with_rinfrc(LLR_2,HD1, max_iter,max_max_iter, mul_mat, add_mat, div_mat,...
                h,str_cn_vn, dc, dev_lsts, nm, v_weights, max_attempt);
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
            BERstat(i0), aver_iter(i0) );
        fprintf(msg)

        msgs  = details_in_lines(ebn0, FER,BER, SER, gen_seq_cnt, K,p, aver_iter, conf_detail, report_fle_nme, save_rslt);
        if save_rslt
            save(report_fle_nme+'.mat');
        end
    end
    
end
