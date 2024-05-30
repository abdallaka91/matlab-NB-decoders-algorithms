clear;%/home/abdallah/Downloads/NB_LDPC_Decoders/call_MV_SF_parforloop.m
comput_SER_BER = false;
ZERO=1; % if  0 then simulate all zeros sequence
plt = 0; % continiously plot FER performance if 1
nm = 4;% V2C m2ssage size
dc1 = [0 1 2];
save_rslt = 1;
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
max_iter = 20;
ebn0 = 3.6; %dB
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
load catched_1000_cap_at_3p6dB_204_102_only.mat
cnt = 0;
h0 = 0*h;
iters_cnt = zeros(1000,1);
parfor jjj = 1 : 1000
    s1 = zeros(1,N);
    y_bin_nse = squeeze(catched_err(:,jjj,:));
    % [info_seq, code_seq, valid_symdrom, y_bin] = generate_and_encode(ZERO, h,G, add_mat, mul_mat, p);
    %             nse = sigma*randn(size(y_bin));
    %             y_bin_nse = y_bin + nse;
    LLRfact = 1;
    unreliable_sat=-inf;
    LLR_20 = 1024*LLR_simple3(y_bin_nse, p,LLRfact , unreliable_sat);
    LLR_20 = round(LLR_20);
    LLR_2 = LLR_20;
    tt=0;
    while tt<3
        s1 = zeros(1,N);
        tt=tt+1;
        [iters, dec_seq, success_dec, synd] = presorted_MVSF_4(LLR_2, max_iter, mul_mat, add_mat, div_mat,...
    h,str_cn_vn, dc, dev_lsts, nm, v_weights);
        iters_cnt(jjj) = iters_cnt(jjj)+iters;
        if success_dec
            cnt = cnt + 1;
            break
        else
            LLR_2 = LLR_20;
            % h1 = h0;
            for i = 1 : M
                if synd(i)
                    % h1(i, str_cn_vn{i}) = 1;
                    s1(str_cn_vn{i}) = 1;
                end
            end
            for j = 1 : N
                if ~s1(j)
                    b = LLR_2(j,:);
                    [~,i1] = max(b);
                    bb=LLR_2(j,dec_seq(j)+1);
                    LLR_2(j,dec_seq(j)+1) = 0;
                    LLR_2(j,i1) = 2*bb;
                end
            end
        end
    end
end
