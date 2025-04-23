8% q=64=> CCSK 0111011001011101011001110000010000101110000111011100100001101011
clear
save_rslt=0;
pth1 = (fullfile(pwd, 'related_functions'));
pth6 = (fullfile(pwd, 'results/'));
addpath(pth1);
rng = 0;
snr_start=-6.5;
snr_end=-6.5;
mainPath = pwd;
filename = fullfile(mainPath, '\bpsk_ccsk_gf64\biawgn_ccsk[6]_gfdim6_n6.plr');
% load(filename)
[p, data] = read_nb_polar_files_parameter(filename);

q=64; CCSK_seq_str = '0111011001011101011001110000010000101110000111011100100001101011';
p=log2(q);
K=42;
UIc = 0;
parforN = 100;
max_gen = 2e5;
max_err_cnt1 = 200; % at low Eb_No(<Eb_No_thrshld)
max_err_cnt2 = 200; %at high Eb_No
SNR_db_thrshld = 3.2;
eta0 = CCSK_seq_str' - 48;
etaq = zeros(q,q);
for i=0:q-1
    etaq(:,i+1)=circshift(eta0,-i);
end
etaqm = 1-2*etaq;
N = length(data(1).channel_sorting);
n=log2(N);
M=N-K;
SNRs_db = arrayfun(@(s) s.SNR, data)';
ii = find(SNRs_db>=snr_start,1);
data(1:ii-1)=[];
SNRs_db = arrayfun(@(s) s.SNR, data)';
ii = find(SNRs_db>snr_end,1);
data(ii:end)=[];
SNRs_db = arrayfun(@(s) s.SNR, data)';
SNRs = 10.^(SNRs_db/10);
N0 = 1./(SNRs);
sigma = sqrt(N0);
I1=bi2de(fliplr(de2bi((0:N-1)', n)));
related_variables_pth = fullfile(mainPath, 'related_variables');
pth3 = (fullfile(related_variables_pth, 'GF_arithm'));
fl_nm = ['arith_' num2str(q) '.mat'];
load(fullfile(pth3, fl_nm));
%%
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
    if SNRs_db(i0)>=SNR_db_thrshld
        max_err_cnt = max_err_cnt2;
    end
    iter_cnt_ = 0;
    FER_ = 0;
    SER_ = 0;
    BER_ = 0;
    gen_seq_cnt_ = 0;
    needed_iters_ = nan(parforN,1);
    msg = sprintf("SNR_dB = %.3f dB, FER = %d/%d = %.8f,// BER = %d/%d = %.8f\n",...
        SNRs_db(i0), FER(i0), gen_seq_cnt(i0), FER(i0)/gen_seq_cnt(i0), BER(i0), gen_seq_cnt(i0)*K*p,...
        BERstat(i0) );
    fprintf(msg)
    sigm =sigma(i0);
    coefs = data(i0).gf_coef;
    coefs=ones(size(coefs));
    chan_idx0 = data(i0).channel_sorting';
    chan_idx=chan_idx0;%(I1+1);
    I=(chan_idx(1:K)); %% unsort if needed
    Ic = (chan_idx(1+K:end)); %% unsort if needed
    info_seq = randi([0 q-1],1,K); 
    u = nan(1, N);
    u(I+1) = info_seq;
    u(Ic+1) = UIc;
    h1= coef_2_coef(coefs,I1);
    x=encode(u, h1,add_mat, mul_mat);
    % u1 = code_decod(x, h1,add_mat, div_mat);
    etax = zeros(q, N);
    for i = 1 : N
        etax(:,i)=etaqm(:,x(i)+1);
    end
    while FER(i0) < max_err_cnt && gen_seq_cnt(i0)<max_gen

        for j = 1 : parforN
            gen_seq_cnt_ = gen_seq_cnt_+1;
            nse = sigm*randn(size(etax));
            y = etax + nse;
            L  = LLR_CCSK(y, etaq, q, N, sigma(i0)^2);
            % [mn,HD_L] = min(L,[],1);
            % HD_L = HD_L-1;

            decw = polar_nb_dec(h1,add_mat, mul_mat, div_mat, L,Ic,UIc, M,N,q,n );
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

        fprintf(repmat('\b',1,length(char(msg))));
        msg = sprintf("SNR_dB = %.3f dB, FER = %d/%d = %.8f,// BER = %d/%d = %.8f\n",...
            SNRs_db(i0), FER(i0), gen_seq_cnt(i0), FER(i0)/gen_seq_cnt(i0), BER(i0), gen_seq_cnt(i0)*K*p,...
            BERstat(i0) );
        fprintf(msg)
        if save_rslt
            save(strcat(report_fle_nme,'.mat'));
        end

    end
end

%%
hold on
Esn0 = SNRs_db;
semilogy(Esn0, FERstat,'-o')
grid on
xlabel("SNR (dB)")
ylabel('FER')
xlim([-20 -8])
ylim([1e-7 1])