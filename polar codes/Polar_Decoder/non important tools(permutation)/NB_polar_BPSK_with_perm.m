% q=64=> CCSK 0111011001011101011001110000010000101110000111011100100001101011
clear
save_rslt=0;
pth1 = (fullfile(pwd, 'related_functions'));
pth6 = (fullfile(pwd, 'results/'));
addpath(pth1);

mainPath = pwd;
filename = fullfile(mainPath, 'data_bpsk\biawgn_gfdim6_n1.mat');
load(filename)

SNRs_db=-5;
N=64;
n=log2(N);
q=64;
p=log2(q);
K=6;
UIc = 0;
parforN = 100;
max_gen = 5e5;
max_err_cnt1 = 200; % at low Eb_No(<Eb_No_thrshld)
max_err_cnt2 = max_err_cnt1; %at high Eb_No
SNR_db_thrshld = 3.2;


% N = length(data(1).channel_sorting);
% n=log2(N);
M=N-K;
% SNRs_db = arrayfun(@(s) s.SNR, data)';
% ii = find(SNRs_db>=snr_start,1);
% data(1:ii-1)=[];
% SNRs_db = arrayfun(@(s) s.SNR, data)';
% ii = find(SNRs_db>snr_end,1);
% data(ii:end)=[];
% SNRs_db = arrayfun(@(s) s.SNR, data)';


SNRs = 10.^(SNRs_db/10);
N0 = 1./(SNRs);
sigma = sqrt(N0);
I1=bi2de(fliplr(de2bi((0:N-1)', n)));
related_variables_pth = fullfile(mainPath, 'related_variables');
pth3 = (fullfile(related_variables_pth, 'GF_arithm'));
fl_nm = ['arith_' num2str(q) '.mat'];
load(fullfile(pth3, fl_nm));
%%
words = (0:q-1);

current_date = datestr(now, 'yyyy_mm_dd');
current_time = datestr(now, 'HH_MM_SS');
[~, report_fle_nme0, ~] = fileparts(filename);
report_fle_nme0 = strcat(report_fle_nme0, '_' ,current_date ,'_' ,current_time);
report_fle_nme = strcat(pth6,report_fle_nme0);

snr_cnt = length(SNRs_db);
FERstat = zeros(snr_cnt,1);
SERstat = zeros(snr_cnt,1);
BERstat = zeros(snr_cnt,1);
gen_seq_cnt = 0*ones(snr_cnt,1);
FER = zeros(snr_cnt,1);
SER = zeros(snr_cnt,1);
BER = zeros(snr_cnt,1);
max_err_cnt = max_err_cnt1;



for i0 = 1 : snr_cnt
    data=load('C:\Users\Abdallah Abdallah\Documents\Personals\matlab-NB-decoders-algorithms\polar codes\Genie_Aided_Coefs_Optim\optimized_coefs\optim_coef_with_perm_-5dB_N64.mat');
    per0 = [56 1 20 55 31 43 2 44 51 29 46 26 37 6 9 48 45 14 27 34 54 24 49 5 0 52 7 41 42 19 28 63 10 36 13 57 17 50 39 30 22 47 32 3 60 8 59 21 35 23 62 16 4 61 40 11 25 58 53 12 15 33 18 38];
    
%     data=load('C:\Users\Abdallah Abdallah\Documents\Personals\matlab-NB-decoders-algorithms\polar codes\Genie_Aided_Coefs_Optim\optimized_coefs\N64_q64_snr-4.99dB_observ500_Entropy_all_param.mat');
%     per0 = mul_mat(1+30,:);

    per0_inv(per0+1)=0:q-1;

    h1=data.h1;

    cnd1=true;
    if SNRs_db(i0)>=SNR_db_thrshld
        max_err_cnt = max_err_cnt2;
    end
    iter_cnt_ = 0;
    FER_ = 0;
    SER_ = 0;
    BER_ = 0;
    gen_seq_cnt_ = 0;
    needed_iters_ = nan(parforN,1);
    msg = sprintf("SNR_dB = %.3f dB,K=%d, FER = %d/%d = %.8f,// BER = %d/%d = %.8f\n",...
        SNRs_db(i0), FER(i0),K, gen_seq_cnt(i0), FER(i0)/gen_seq_cnt(i0), BER(i0), gen_seq_cnt(i0)*K*p,...
        BERstat(i0) );
    %     fprintf(msg)
    sigm =sigma(i0);
    %     coefs = data(i0).gf_coef;
    chan_idx0 = data.ch_idx'-1;
    chan_idx=chan_idx0;%(I1+1);
    I=sort(chan_idx(1:K)); %% unsort if needed
    Ic = sort(chan_idx(1+K:end)); %% unsort if needed
    info_seq = randi([0 q-1],1,K);
    u = nan(1, N);
    u(I+1) = info_seq;
    u(Ic+1) = UIc;
    %     h1= coef_2_coef(coefs,I1); %
    %      data1=load('C:\Users\Abdallah Abdallah\Documents\Personals\matlab-NB-decoders-algorithms\polar codes\Genie_Aided_Coefs_Optim\optimized_coefs\N64_q64_snr0.01dB_observ1000_Entropy.mat');
    x=encode_coef_perm(u, h1,add_mat, mul_mat, per0);
    u1 = code_decod_perm(x, h1,add_mat, div_mat, per0);
    alph_bin =  fliplr(dec2bin(words, p) - 48);
    alph_bin_mod = (-1).^alph_bin;
    etax = alph_bin_mod(x+1, :);
    while FER(i0) < max_err_cnt && gen_seq_cnt(i0)<max_gen
        parfor j = 1 : parforN
            gen_seq_cnt_ = gen_seq_cnt_+1;
            nse = sigm*randn(size(etax));
            y = etax + nse;
            L = -LLR_BPSK_GFq_2D(y, sigm, alph_bin);
            [mn,HD_L] = min(L,[],1);
            HD_L = HD_L-1;
            decw = polar_nb_dec_with_perm(h1,add_mat, mul_mat, div_mat, per0,per0_inv, L,Ic,UIc, M,N,q,n );
            rec_info_seq = decw(I+1);
            nd = sum(rec_info_seq~=info_seq);
            if  nd ~=0
                FER_ = FER_ +1;
                [SER0, BER0] = SER_BER(info_seq,rec_info_seq,p, 0, 0);
                SER_ = SER0 + SER_;
                BER_ = BER0 + BER_;
            end
        end
        SER(i0) = SER_;
        BER(i0) = BER_;
        gen_seq_cnt(i0) = gen_seq_cnt_;

        FER(i0) = FER_;
        FERstat(i0)=FER(i0)/gen_seq_cnt(i0);
        SERstat(i0)=SER(i0)/(gen_seq_cnt(i0)*K);
        BERstat(i0)=BER(i0)/(gen_seq_cnt(i0)*K*p);

        if ~cnd1
            fprintf(repmat('\b',1,length(char(msg))));
        else
            cnd1=false;
        end
        msg = sprintf("SNR_dB = %.3f dB, K=%d, FER = %d/%d = %.8f,// BER = %d/%d = %.8f\n",...
            SNRs_db(i0),K,  FER(i0), gen_seq_cnt(i0), FER(i0)/gen_seq_cnt(i0), BER(i0), gen_seq_cnt(i0)*K*p,...
            BERstat(i0) );
        fprintf(msg)
        if save_rslt
            save(strcat(report_fle_nme,'.mat'));
        end

    end
end

%%
% Esn0 = SNRs_db;
% semilogy(Esn0, FERstat,'-o')
% grid on
% xlabel("SNR (dB)")
% ylabel('FER')
% xlim([-1.5 1])
% ylim([1e-4 1])