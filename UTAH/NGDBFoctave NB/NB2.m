clear

ZERO = 0; % all zeros seq
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
K = N-M;
Nb = N*p;
Kb = K*p;
% K=175;
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
Ebnos           = 4.4;  % SNR values to simulate
snr_cnt         = length(Ebnos);
max_gen           = 1e6;
F                 = 100;      % Number of frame errors to observe
Fchan             = 1;    % Number of channel frames to generate per batch
R                 = 0.5;  % code rate
T                 = 1000;    % Max iterations for decoding
w                 = 1;   % Syndrome weight parameter
Ymax              = 295;    % Channel sample saturation magnitude
theta             = -2.2;   % Flipping threshold
eta               = 1.0;     % Perturbation noise scale parameter


GDBFerrors       = zeros(snr_cnt, 1);
GDBFframeerrors  = zeros(snr_cnt, 1);
GDBFframeundetected  = zeros(snr_cnt, 1);
NFrames          = zeros(snr_cnt, 1);

p1 = 1;
Rate = p1*K/N; %p1 is nb of bits per channel use with the modulation, for example for bpsk it is 1
ebn0_n = 10.^(Ebnos/10);
N0 = 1./(Rate*ebn0_n);
sigma = sqrt(N0/2);

snr = -10*log10(2*sigma.^2);
%%
rng(1); % repetitive noise generation
alph_bin =  fliplr(dec2bin(words, p) - 48);
alph_bin_mod = (-1).^alph_bin;

NFrames(snr_cnt)=0;
[G,~] = Generator_matrix_G_from_full_rank_H(h, add_mat, mul_mat, div_mat);

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
iter_cnt = zeros(snr_cnt,1);

for sdx=1:snr_cnt
    csigma = sigma(sdx);
    nsigma    = eta*csigma;
    msg = sprintf("EbNo = %.3f dB, FER = %d/%d = %.8f,// BER = %d/%d = %.8f, aver_iter = %.3f\n",...
            Ebnos(sdx), FER_cnt(sdx), F_gen(sdx), FERR(sdx), BER_cnt(sdx),B_gen(sdx),...
            BERR(sdx), aver_iter(sdx) );

        fprintf(msg)
        keepGoing=1;
    while (keepGoing)
        

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



        noisevec = csigma*randn(Fchan,Nb);              % Matrix of AWGN noise
        y        = min(max(y_bin_mod1D + noisevec,-Ymax),Ymax);  % Clipped channel samples
        dB = (1-sign(y))/2;
        nerrb_decs = sum(dB~=y_bin1D);
        ber_decs = nerrb_decs/Nb;
        dgf = HD_bin_gf(sign(y),p);
        nerrS_decs = sum(dgf~=code_seq);
        Ser_decs = nerrS_decs/N;
        [d, failed, Sgf, Sb, E, seqb, seqgf, iters] = decodeGDBFvecNB2(p, str_cn_vn,mul_mat, add_mat, y,hl, h, N, T, w, theta, nsigma);
        rec_info_seq = seqgf(1:K);
        
        if failed~=0
        [nerrS_decd, nerrB_decd] = SER_BER(info_seq,rec_info_seq,p, 0, 0);
        else
            nerrS_decd=0;
            nerrB_decd=0;
        end


        ber_decd = nerrB_decd/Kb;
        Ser_decd = nerrS_decd/K;

        nerrF_decd=0;
        if Ser_decd>0
            nerrF_decd=1;
        end

        SER_cnt(sdx) = SER_cnt(sdx)+nerrS_decd;
        BER_cnt(sdx) = BER_cnt(sdx)+nerrB_decd;
        FER_cnt(sdx) = FER_cnt(sdx)+nerrF_decd;

        S_gen(sdx) = S_gen(sdx)+K;
        B_gen(sdx) = B_gen(sdx)+Kb;
        F_gen(sdx) = F_gen(sdx)+1;

        SERR(sdx) = SER_cnt(sdx)/S_gen(sdx) ;
        BERR(sdx) = BER_cnt(sdx)/B_gen(sdx) ;
        FERR(sdx) = FER_cnt(sdx)/F_gen(sdx) ;

        iter_cnt(sdx) = iter_cnt(sdx) + iters;
        aver_iter(sdx) = iter_cnt(sdx)/F_gen(sdx);

        fprintf(repmat('\b',1,length(char(msg))));
        msg = sprintf("EbNo = %.3f dB, FER = %d/%d = %.8f,// BER = %d/%d = %.8f, aver_iter = %.3f\n",...
            Ebnos(sdx), FER_cnt(sdx), F_gen(sdx), FERR(sdx), BER_cnt(sdx),B_gen(sdx),...
            BERR(sdx), aver_iter(sdx) );
        fprintf(msg)

        if (FER_cnt(sdx) > F || F_gen(sdx)>max_gen)
            keepGoing = 0;
        end

        

    end
end