clear;%/home/abdallah/Downloads/NB_LDPC_Decoders/call_MV_SF_parforloop.m
restoredefaultpath
comput_SER_BER = false;

% v0v1 = [2.4 1];
v0v1 = [0.5 0.25];
vi = 3;
nui = 2; % nb of locations to change
L=6;
LLRfact = 4;
unreliable_sat=-50;

max_err_cnt = 50;
max_gen = 0.5e6;
max_iter = 10;
ebn0 = 4.3; %dB


pth1 = (fullfile(pwd, 'related_functions'));
addpath(pth1);
pth2 = (fullfile(pwd, 'related_variables'));
pth3 = (fullfile(pwd, 'related_variables/GF_arithm'));
pth4 = (fullfile(pwd, 'related_variables/alists'));
pth5 = (fullfile(pwd, 'related_variables/alists/matrices'));
pth6 = (fullfile(pwd, 'results/'));
pth7 = (fullfile(pwd, 'helpful functions/'));
addpath(pth7);


% H_matrix_mat_fl_nm = '837_124_32';
% H_matrix_mat_fl_nm = '204.102.3.6.16';
H_matrix_mat_fl_nm = '837_124_32';
load([fullfile(pth4, H_matrix_mat_fl_nm) '.mat']);
H=h;
h=H; %%--------------------->>>>>> K is not here size(H,2)-size(H,1)
K = 726;
N = size(h,2);
M = size(h,1);

p = 5;%ceil(log2(max(max(H))+0.1));
q = 2^p;
words = (0:q-1);

fl_nm = ['arith_' num2str(q) '.mat'];
if  exist(fullfile(pth3, fl_nm), 'file') == 2
    load(fullfile(pth3, fl_nm));
else
    add_mat = GF_arithm_matrix(q, 'add');
    mul_mat = GF_arithm_matrix(q, 'mul');
    div_mat = GF_arithm_matrix(q, 'div');
    save(fullfile(pth3, ['arith_' num2str(q) '.mat']), 'add_mat' ,'mul_mat','div_mat')
end

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
info_seq = 0*randi([0 q-1], 1, K);
% info_seq_bit=(fliplr(de2bi(info_seq,p)))';
info_seq_bit = fliplr(dec2bin(info_seq, p) - 48)';
info_seq_bit=info_seq_bit(:);
code_seq = zeros(1,N);
% y_bin0 = de2bi(code_seq,p);
y_bin0 = fliplr(dec2bin(code_seq, p) - 48);
y_bin = (-1).^y_bin0;
%%
alph1 = 1:vi;
combs = alph1;
for j0 = 2:nui
    % combs = combvec(combs, alph1);
    combs = CombVec(combs, alph1);
end
combs = combs';

snr_cnt = length(sigma);
FERstat = zeros(snr_cnt,1);
gen_seq_cnt = zeros(snr_cnt,1);
FER = zeros(snr_cnt,1);


parforN = 100;

for i0 = 1 : snr_cnt
    FER_ = 0;
    gen_seq_cnt_ = 0;
    msg = sprintf("EbNo = %.3f dB, V0_V1= [%.3f %.3f], Error frames/Total frames = %d/%d => FER = %.8f\n",...
        ebn0(i0), v0v1(1), v0v1(2), FER(i0), gen_seq_cnt(i0), FER(i0)/gen_seq_cnt(i0));
    fprintf(msg)
    sigm =sigma(i0);
    while FER(i0) < max_err_cnt && gen_seq_cnt(i0)<max_gen
        tic
        parfor j = 1 : parforN
            gen_seq_cnt_ = gen_seq_cnt_+1;
            nse = sigm*randn(size(y_bin));
            y_bin_nse = y_bin + nse;
%             LLR_2 = LLR_BPSK_GFq_2D(y_bin_nse, sigm);
            LLR_2 = LLR_simple3(y_bin_nse, p,LLRfact , unreliable_sat);
            [iter, dec_seq, success_dec] = nb_ldpc_MV_SF_opt1(LLR_2, vi,v0v1, nui,L, max_iter, mul_mat, add_mat, div_mat,combs, h);
            rec_info_seq = dec_seq(N+1-K:N);
            nd = sum(rec_info_seq~=info_seq);
            if nd ~=0
                FER_ = FER_ +1;
            end
        end
        TOC=toc;
        gen_seq_cnt(i0) = gen_seq_cnt_;
        FER(i0) = FER_;
        FERstat(i0)=FER(i0)/gen_seq_cnt(i0);

        fprintf(repmat('\b',1,length(char(msg))));
        msg = sprintf("EbNo = %.3f dB, V0_V1= [%.3f %.3f], Error frames/Total frames = %d/%d => FER = %.8f\n",...
            ebn0(i0), v0v1(1), v0v1(2), FER(i0), gen_seq_cnt(i0), FER(i0)/gen_seq_cnt(i0));
        fprintf(msg)

    end
end
%%
% ebn00 = ebn0';
% FERstat0 = FER./gen_seq_cnt;
% 
% FERstat1 = [0.09 0.04 0.015 0.005 0.0015 0.00042 0.0 0.0];
% ebn01 = 4:0.1:4.7;
% 
% figure(1)
% 
% 
% semilogy(ebn00, FERstat0,'bo:', 'LineWidth',1.2)
% hold on
% semilogy(ebn01, FERstat1,'ro:', 'LineWidth',1.2)
% xlabel('E_b/N_0 (dB)')
% ylabel('FER (Log scale)')
% grid on
% xlim([3 5])
% ylim([10e-5 1])
% title('GF(256), (255, 175), L=4')
% legend({'our code, itr=10, v0=0.5,v1=0.25'; 'author code, itr=15, v0v1v2'})
%%

% ebn00 = ebn0';
% FERstat0 = FER./gen_seq_cnt;
% 
% FERstat0 = [91/400 94/800 80/2000 80/5800 80/23700 95/165500 8/83000];
% ebn00 = 4:0.1:4.6;
% 
% FERstat1 = [0.2 0.1 0.03 0.007 0.0015 0.00025 0.000055];
% ebn01 = 4:0.1:4.6;
% 
% figure(1)
% 
% 
% semilogy(ebn00, FERstat0,'bo:', 'LineWidth',1.2)
% hold on
% semilogy(ebn01, FERstat1,'ro:', 'LineWidth',1.2)
% xlabel('E_b/N_0 (dB)')
% ylabel('FER (Log scale)')
% grid on
% xlim([3 5])
% ylim([10e-6 1])
% title('GF(32), (837, 726), L=4')
% legend({'our code, L=4, itr=10, v0=0.5,v1=0.25'; 'author code, L=4, itr=10, v0=0.5 v1=0.25'})
%%
% FERstat1 = [0.3 0.21 0.12 0.055 0.025 0.008 0.003 0.0008 0.00025 0.000065];
% ebn01 = 4:0.1:4.9;
% 
% figure(3)
% hold on
% semilogy(ebn01, FERstat1,'ro:', 'LineWidth',1.2)
% xlabel('E_b/N_0 (dB)')
% ylabel('FER (Log scale)')
% grid on
% xlim([3 5])
% ylim([10e-5 1])
% title('GF(256), (255, 175), L=4')
% legend({'our code, itr=10, v=1'; 'author code, itr=15, v=1'})




