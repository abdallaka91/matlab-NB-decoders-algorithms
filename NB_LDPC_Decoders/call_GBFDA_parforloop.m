clear

pth1 = (fullfile(pwd, 'related_functions\'));
addpath(pth1);
pth2 = (fullfile(pwd, 'related_variables\'));
pth3 = (fullfile(pwd, 'related_variables\GF_arithm'));
pth4 = (fullfile(pwd, 'related_variables\alists\'));
pth5 = (fullfile(pwd, 'related_variables\alists\matrices\'));
pth6 = (fullfile(pwd, 'results\'));

H_matrix_mat_fl_nm = 'EG_255_175';
load([fullfile(pth4, H_matrix_mat_fl_nm) '.mat']); h=H;
H = full(h);
% h=H; %%--------------------->>>>>> K is not here size(H,2)-size(H,1)
K = 175;
N = size(h,2);
M = size(h,1);

ebn0 = 4.0:0.1:4.9;
p = 8;%ceil(log2(max(max(H))+0.1));%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q = 2^p;
max_err_cnt = 50;
max_gen = 1e6;
v = 1;
max_iter = 10;
max_iter_needed = 1;
LLRfact = 0.5;


fl_nm = ['arith_' num2str(q) '.mat'];
if  exist(fullfile(pth3, fl_nm), 'file') == 2
    load(fullfile(pth3, fl_nm));
else
    add_mat = GF_arithm_matrix(q, 'add');
    mul_mat = GF_arithm_matrix(q, 'mul');
    div_mat = GF_arithm_matrix(q, 'div');
    save(fullfile(pth3, ['arith_' num2str(q) '.mat']), 'add_mat' ,'mul_mat','div_mat')
end


lgnd = ['GBFDA__' H_matrix_mat_fl_nm '_' num2str(max_iter) '_Iter_V_'...
    num2str(v)];

% N = size(H,2);
% M = size(H,1);
% K = N-M;
p1 = 1;
Rate = p1*K/N; %p1 is nb of bits per channel use with the modulation, for example for bpsk it is 1

snr_1 = 10.^(ebn0/10);
N00 = 1./snr_1;
N0 = N00/Rate;
sigma = sqrt(N0/2);


%%
info_seq = 0*randi([0 q-1], 1, K);
info_seq_bit=(fliplr(de2bi(info_seq,p)))';
info_seq_bit=info_seq_bit(:);
code_seq = zeros(1,N);
y_bin0 = de2bi(code_seq,p);
y_bin = (-1).^y_bin0;
H = sparse(H);

%%
Wmn = cell(M,1);
syndrm = zeros(1,M);

lst1 = cell(M,1);
dc = zeros(M,1);
for i  =1 : M
    lst1{i} = find(H(i,:));
    dc(i) = length(lst1{i});
end

snr_cnt = length(sigma);
FERstat = zeros(snr_cnt,1);
gen_seq_cnt = zeros(snr_cnt,1);
FER = zeros(snr_cnt,1);


parforN = 100;

for i0 = 1 : snr_cnt
    FER_ = 0;
    gen_seq_cnt_ = 0;
    msg = sprintf("EbNo = %.3f dB, V = %.3f, Error frames/Total frames = %d/%d => FER = %.8f\n",...
        ebn0(i0), v, FER(i0), gen_seq_cnt(i0), FER(i0)/gen_seq_cnt(i0));
    fprintf(msg)
    sigm =sigma(i0);
    while FER(i0) < max_err_cnt && gen_seq_cnt(i0)<max_gen
        tic
        parfor j = 1 : parforN
            gen_seq_cnt_ = gen_seq_cnt_+1;
            nse = sigm*randn(size(y_bin));
            y_bin_nse = y_bin + nse;

%             LLR_2 = LLR_BPSK_GFq_2D(y_bin_nse, sigm);

            LLR_2 = LLRfact*LLR_simple2(y_bin_nse,p, sigm);

            [needed_iter, dec_seq, is_code] = Serial_Enhanced_GBFDA(LLR_2, dc, lst1, q, h,max_iter, add_mat, mul_mat, div_mat, v);
            %             rec_info_seq = dec_seq(N+1-K:N);
            %             nd = sum(rec_info_seq~=info_seq);
            nd = sum(code_seq~=dec_seq);
            if nd ~=0
                FER_ = FER_ +1;
            end
        end
        TOC=toc;
        gen_seq_cnt(i0) = gen_seq_cnt_;
        FER(i0) = FER_;
        FERstat(i0)=FER(i0)/gen_seq_cnt(i0);

        fprintf(repmat('\b',1,length(char(msg))));
        msg = sprintf("EbNo = %.3f dB, V = %.3f, Error frames/Total frames = %d/%d => FER = %.8f\n",...
            ebn0(i0), v, FER(i0), gen_seq_cnt(i0), FER(i0)/gen_seq_cnt(i0));
        fprintf(msg)

    end
end

%%
figure(1)

FERstat0 = FER./gen_seq_cnt;
semilogy(ebn0, FERstat0,'bo:', 'LineWidth',1.2)
hold on
xlabel('E_b/N_0 (dB)')
ylabel('FER (Log scale)')
grid on
xlim([4 5.5])
ylim([10e-6 1])
legend('(255,175)NB-LDPC, MV-SF, it=10, GF(256)')



