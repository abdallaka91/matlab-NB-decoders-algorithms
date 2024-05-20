clear

pth1 = (fullfile(pwd, 'related_functions\'));
addpath(pth1);
pth2 = (fullfile(pwd, 'related_variables\'));
pth3 = (fullfile(pwd, 'related_variables\GF_arithm'));
pth4 = (fullfile(pwd, 'related_variables\alists\'));
pth5 = (fullfile(pwd, 'related_variables\alists\matrices\'));
pth6 = (fullfile(pwd, 'results\'));

H_matrix_mat_fl_nm = 'generated_6x12_GF4';

load([fullfile(pth4, H_matrix_mat_fl_nm) '.mat']);
H = h;

refresh_figure_every = 20;
comput_SER_BER = true;
ebn0 = 0:0.2:5.20;
p = ceil(log2(max(max(H))+0.1));
q = 2^p;
max_err_cnt = 50;
max_gen_seq = 3e4;
v = 1;
max_iter = 50;
max_iter_needed = 1;



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

N = size(H,2);
M = size(H,1);
K = N-M;
Rate = K/N;


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
h_gf = gf(full(H),p);
%%
syndrm = zeros(1,M);

cnt = 0;
lst1 = cell(M,1);
dc = zeros(M,1);
for i  =1 : M
    lst1{i} = find(H(i,:));
    dc(i) = length(lst1{i});
end

lst_comb_all = cell(M,1);
for i = 1 : M
    temp_gf2 = gf(zeros(1, N), p);
    idx1 = lst1{i,1};
    l1 = length(idx1)-1;
    lst_comb = nan(q^(l1), length(idx1));
    %     for j1 = 1 : length(idx1)
    j1 = 1;
    idx2 = idx1;
    idx2(j1) = [];

    temp = zeros(1, l1);
    temp(1) = -1;
    for j2 = 1 : q^l1
        temp = comb_q_m(temp, q, l1);
        temp_gf = gf(temp,p);
        arr_coef = h_gf(i,idx2);
        temp_gf1 = temp_gf.*arr_coef;
        compl = gf(0,p);
        for j3 = 1 : l1
            compl = compl+temp_gf1(j3);
            cnt = cnt+1;
        end
        compl1 = compl/h_gf(i,idx1(j1));

        temp_gf2(idx2) = temp_gf;
        temp_gf2(idx1(j1)) = compl1;
        valid = temp_gf2*h_gf(i,:)';
        lst_comb(j2,:) = ([double(compl1.x) temp]);
    end
    lst_comb_all{i,1 } = lst_comb;
end

% i10 = 5;dd=gf(lst_comb_all{i10,1 },p)*h_gf(i10,lst1{i10,1})'; 

%%
snr_cnt = length(sigma);
BERstat = zeros(snr_cnt,1);
SERstat = zeros(snr_cnt,1);
FERstat = zeros(snr_cnt,1);
gen_sym_cnt = zeros(snr_cnt,1);
gen_bit_cnt = zeros(snr_cnt,1);
gen_seq_cnt = zeros(snr_cnt,1);
FER = zeros(snr_cnt,1);
SER = zeros(snr_cnt,1);
BER = zeros(snr_cnt,1);
FER_HD = zeros(snr_cnt,1);
SER_HD = zeros(snr_cnt,1);
BER_HD = zeros(snr_cnt,1);
BER_HDstat = zeros(snr_cnt,1);
SER_HDstat = zeros(snr_cnt,1);
FER_HDstat = zeros(snr_cnt,1);

for i0 = 1 : snr_cnt
    last_refresh_cnt = 1;
    while FER(i0)<max_err_cnt && gen_seq_cnt(i0)<max_gen_seq
        gen_seq_cnt(i0) = gen_seq_cnt(i0)+1;
        gen_bit_cnt(i0) = gen_bit_cnt(i0) + K*p;
        gen_sym_cnt(i0) = gen_sym_cnt(i0) + K;


        nse = sigma(i0)*randn(size(y_bin));
        y_bin_nse = y_bin + nse;

        %         bin_HD1 = ones(size(y_bin_nse));
        %         bin_HD1(y_bin_nse>=0)=0;

        LLR_2 = LLR_BPSK_GFq_2D(y_bin_nse, sigma(i0));

        [needed_iter, dec_seq, is_code] = NB_min_sum(max_iter, LLR_2, dc, lst_comb_all, lst1, q, h, h_gf);

        if needed_iter>max_iter_needed
            max_iter_needed = needed_iter;
        end

        rec_info_seq = dec_seq(N+1-K:N);
        nd = sum(rec_info_seq~=info_seq);
        if nd ~=0
            FER(i0) = FER(i0)+1;
            SER(i0)=SER(i0)+nd;
        end

        FERstat(i0)=FER(i0)/gen_seq_cnt(i0);
        SERstat(i0)= SER(i0)/gen_sym_cnt(i0);

        if comput_SER_BER

            rec_info_seq_bit=(fliplr(de2bi(rec_info_seq,p)))';
            rec_info_seq_bit=rec_info_seq_bit(:);

            BER(i0)=BER(i0)+size(find(info_seq_bit~=rec_info_seq_bit),1);
            BERstat(i0)=BER(i0)/gen_bit_cnt(i0);

            [~,rec_seq_HD]=max(LLR_2, [],1);
            rec_seq_HD = rec_seq_HD-1;
            rec_info_HD = rec_seq_HD(N+1-K:N);
            nd = sum(rec_info_HD~=info_seq);

            if nd ~=0
                FER_HD(i0) = FER_HD(i0)+1;
                SER_HD(i0)=SER_HD(i0)+nd;
                FER_HDstat(i0)=FER_HD(i0)/gen_seq_cnt(i0);
                SER_HDstat(i0)= SER_HD(i0)/gen_sym_cnt(i0);
            end
            rec_info_seq_bit_HD=(fliplr(de2bi(rec_info_HD,p)))';
            rec_info_seq_bit_HD=rec_info_seq_bit_HD(:);
            BER_HD(i0)=BER_HD(i0)+size(find(info_seq_bit~=rec_info_seq_bit_HD),1);
            BER_HDstat(i0)=BER_HD(i0)/gen_bit_cnt(i0);
        end

        if gen_seq_cnt(i0)==last_refresh_cnt

            last_refresh_cnt = last_refresh_cnt+refresh_figure_every;
            figure(1)

            semilogy(ebn0, FERstat,'ro:', 'LineWidth',1.2)
            hold on
            xlabel('E_b/N_0 (dB)')
            ylabel('FER (Log scale)')
            grid on
            %             legend(lgnd, 'Interpreter','none')
            title(lgnd, 'Interpreter','none')

            hold off
            xlim([ebn0(1) ebn0(end)+1])
            ylim([1e-6 2])
            pause(0.2)
        end
    end
    Whos = whos;
    workspaceInfo = Whos;
    workspaceStruct = struct();
    for hh = 1:length(workspaceInfo)
        varName = workspaceInfo(hh).name;
        workspaceStruct.(varName) = eval(varName);
    end
    save(fullfile(pth6,[lgnd '.mat']), 'workspaceStruct')
    saveas(gcf, fullfile(pth6,[lgnd '.fig']))
end
%%

