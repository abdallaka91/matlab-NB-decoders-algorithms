clear;%/home/abdallah/Downloads/NB_LDPC_Decoders/call_MV_SF_parforloop.m
comput_SER_BER = false;
ZERO=0;
plt = 0;
nm = 4;
nc = 3*nm;
c2v_comp_fact=50;

max_err_cnt = 100;
max_gen = 1e6;
max_iter = 10;
ebn0 = 4; %dB
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
end


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
parforN = 100;
for i0 = 1 : snr_cnt
    iter_cnt_ = 0;
    FER_ = 0;
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
            LLR_2 = -1024*LLR_simple3(y_bin_nse, p,LLRfact , unreliable_sat);
            LLR_2 = round(LLR_2);
            [iter, dec_seq, success_dec] = EMS1(...
                LLR_2', max_iter, mul_mat, add_mat, div_mat, h,str_cn_vn, dc,str_vn_cn, dv, nm, nc, c2v_comp_fact);
            iter_cnt_ = iter_cnt_ + iter;
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
        
    end
end



%%
clear conf_detail
conf_detail.a1fl_nme = sprintf("H matrix : %s",H_matrix_mat_fl_nm);
conf_detail.a2Code = sprintf("N = %d, M = %d, K = %d, GF(%d)",N,M,K,q);
conf_detail.a3votesW = sprintf('V0 = %f, V1 = %f', v_weights(1) , v_weights(2));
conf_detail.a4iter = sprintf("H matrix : %d",max_iter);
conf_detail.a5max_seq = sprintf("max seq generation : %d", max_gen);
conf_detail.a6fl_nme = sprintf("max error error frame detection: %d",max_err_cnt);
conf_detail.a7dev_cnt = sprintf("nb of defiation paths : %d", length(dev_lsts_i));
conf_detail.a8reg_widths = sprintf("regions width (high to low reliable) : {") + sprintf(" %d ", dc11) + sprintf("}");
formatted_str = '';
for i = 1:length(di)
    formatted_str = [formatted_str, sprintf('[%d %d]', di{i}(2:end))];
    if i < length(di)
        formatted_str = [formatted_str, ' || '];
    end
end
conf_detail.a9di1_di2 = sprintf("%s\n", formatted_str);


msg  = details_in_lines(ebn0, FER,BER, SER, gen_seq_cnt, K,p, aver_iter, conf_detail);
msg{end}


%%
function plot_FER(b, ebn01, H_matrix_mat_fl_nm)
a=[40./[48 55 81 108 171 440 867 2007 6238 13415 25430 60105 149700 403359 1071371 02352556 04296652 11037931] 35/23048410];
ebn00 = 1.4:0.2:5;
figure(1)
% semilogy(ebn00, a,'bo:', 'LineWidth',1.2
semilogy(ebn01, b,'ro:', 'LineWidth',1.2)
xlabel('E_b/N_0 (dB)')
ylabel('BER (Log scale)')
grid on
xlim([0 6])
ylim([10e-7 1])
title(H_matrix_mat_fl_nm)
% legend({'EMS (Cedric code) it=15, n-vc=n-cv=4'; ['Presorted-Syndrome-based-Request' newline  'Modified MVSF, iter=15, deviat=9']})
pause(0.01)
end
%%

