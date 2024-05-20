clear;%/home/abdallah/Downloads/NB_LDPC_Decoders/call_MV_SF_parforloop.m
comput_SER_BER = false;
ZERO = 1;
plt = 0;
nm = 4;
dc1 = [0 1 2];
Dev_pos_cnt = 4;%length(dc1)-1;
di = cell(length(dc1),1);
di{1} = [0 0];
di{2} = [0 1];
di{3} = [0 2 1];


% dc1 = [0 2];
% di = cell(length(dc1),1);
% di{1} = [0 0];
% di{2} = [0 2 1];


max_err_cnt = 100;
max_gen = 1e6;
max_iter = 15;
ebn0 = 3; %dB
p = 6;
q = 2^p;
K = 8;

pth1 = (fullfile(pwd, 'related_functions'));
addpath(pth1);
pth2 = (fullfile(pwd, 'related_variables'));
pth3 = (fullfile(pwd, 'related_variables/GF_arithm'));
pth4 = (fullfile(pwd, 'related_variables/alists'));
pth5 = (fullfile(pwd, 'related_variables/alists/matrices'));
pth6 = (fullfile(pwd, 'results/'));


words = (0:q-1);

H_matrix_mat_fl_nm = 'N96_K48_GF64_non_exponen_form';
load([fullfile(pth4, H_matrix_mat_fl_nm) '.mat']);
% h=H;
h = full(h);

N = size(h,2);
M = size(h,1);

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

% M = size(H, 1);
% N = size(H, 2);
% K = N-M;

p1 = 1;
Rate = p1*K/N; %p1 is nb of bits per channel use with the modulation, for example for bpsk it is 1

ebn0_n = 10.^(ebn0/10);
N0 = 1./(Rate*ebn0_n);
sigma = sqrt(N0/2);
snr = -10*log10(2*sigma.^2);

%%
% info_seq = zeros(1,K);
% info_seq_bit = fliplr(dec2bin(info_seq, p) - 48)';
% info_seq_bit=info_seq_bit(:);
% code_seq = zeros(1, N);
% y_bin0 = fliplr(dec2bin(code_seq, p) - 48);

alph_bin =  fliplr(dec2bin(words, p) - 48);
alph_bin_mod = (-1).^alph_bin;

[G,~] = Generator_matrix_G_from_full_rank_H(h, add_mat, mul_mat, div_mat);
[info_seq, code_seq, valid_symdrom, y_bin] = generate_and_encode(1, h,G, add_mat, mul_mat, p);
%%
% vi = 3;
% nui=2;
% L = vi^nui;
% alph1 = 1:vi;
% combs = alph1;
% for j0 = 2:nui
%     % combs = combvec(combs, alph1);
%     combs = CombVec(combs, alph1);
% end
% combs = combs';


str_vn_cn = cell(N,1);
dv = zeros(N,1);
for j = 1 : N
    str_vn_cn{j, 1} = (find(h(:,j)))';
    dv(j) = length(str_vn_cn{j});
end

tmp0 = 160;
% tmp1 = 0:0.2:2;
v_weights_m = 0:4:1024;

v_weights_m = CombVec(v_weights_m,tmp0);
% v_weights_m = CombVec(v_weights_m,tmp1);
v_weights_m = v_weights_m';

snr_cnt = length(v_weights_m);
FERstat = zeros(snr_cnt,1);
aver_iter = zeros(snr_cnt,1);
gen_seq_cnt = eps*ones(snr_cnt,1);
FER = zeros(snr_cnt,1);
iter_cnt = zeros(snr_cnt,1);


parforN = 100;

v_weights = [0.4 0.4 0.18];
for i0 = 1 : snr_cnt
    v_weights = v_weights_m(i0,:);
    % v_weights = [400 200 200];
    iter_cnt_ = 0;
    FER_ = 0;
    gen_seq_cnt_ = 0;
    msg = sprintf("%d: %.3f %.3f, Error frames/Total frames = %d/%d => FER = %.8f, aver_iter = %.3f\n",...
        i0, v_weights(1),v_weights(2), FER(i0), gen_seq_cnt(i0), FER(i0)/gen_seq_cnt(i0), 0);
    fprintf(msg)
    sigm =sigma;
    while FER(i0) < max_err_cnt && gen_seq_cnt(i0)<max_gen
        tic
        parfor j = 1 : parforN
            [info_seq, code_seq, valid_symdrom, y_bin] = generate_and_encode(ZERO, h,G, add_mat, mul_mat, p);
            gen_seq_cnt_ = gen_seq_cnt_+1;
            nse = sigm*randn(size(y_bin));
            y_bin_nse = y_bin + nse;
            LLRfact = 1;
            unreliable_sat=-inf;
            LLR_2 = 1024*LLR_simple3(y_bin_nse, p,LLRfact , unreliable_sat);
            LLR_2 = round(LLR_2);
            [iter, dec_seq, success_dec] = presorted_MVSF_2(LLR_2', max_iter, mul_mat, add_mat, div_mat, h,str_cn_vn, dc, ...
                dev_lsts,dev_pos,nm,v_weights, Dev_pos_cnt);
            iter_cnt_ = iter_cnt_ + iter;
            rec_info_seq = dec_seq(1:K);
            nd = sum(rec_info_seq~=info_seq);
            if  nd ~=0
                FER_ = FER_ +1;
            end
        end
        gen_seq_cnt(i0) = gen_seq_cnt_;
        iter_cnt(i0) = iter_cnt_;
        aver_iter(i0) = iter_cnt(i0)/gen_seq_cnt(i0);
        TOC=toc;
        
        FER(i0) = FER_;
        FERstat(i0)=FER(i0)/gen_seq_cnt(i0);

        fprintf(repmat('\b',1,length(char(msg))));
        msg = sprintf("%d: %.3f %.3f, Error frames/Total frames = %d/%d => FER = %.8f, aver_iter = %.3f\n",...
        i0, v_weights(1),v_weights(2), FER(i0), gen_seq_cnt(i0), FER(i0)/gen_seq_cnt(i0), aver_iter(i0));
        fprintf(msg)

    end
end

%%

msg='';
for tt=1:snr_cnt
    formatSpec = msg + "%s";
    msg0 =  sprintf("EbNo = %.3f dB, Error frames/Total frames = %d/%d => FER = %.8f, aver_iter = %.3f\n",...
            ebn0(tt), FER(tt), gen_seq_cnt(tt), FER(tt)/gen_seq_cnt(tt), iter_cnt(tt)/gen_seq_cnt(tt));
    msg = sprintf(formatSpec,msg0);
end
fprintf(msg)
