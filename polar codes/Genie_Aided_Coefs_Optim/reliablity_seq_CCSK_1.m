clear
% rng(1);
pth1 = (fullfile(pwd, 'related_functions'));
addpath(pth1);
mainPath = pwd;
q=256;
SNRs_db = -10:0.5:-10;
N = 64;
n_obs=500;
low2highreliab = 1;

if q==64
    CCSK_seq_str = '0111011001011101011001110000010000101110000111011100100001101011'; %GF64
elseif q==128
    CCSK_seq_str = '10001010101000001101100110101101000001100011011001101100000101000011011110011100101010000111101001001110111111111000110111010100';%GF128
elseif q==256
    CCSK_seq_str = '0111010101010010010011011010111100100110001000010100100100011010101000010101011101101110000001100110011100001100100101100011100010011110011100101011110100111001100010110000101100000100100111011110000110000111010100010111111101111111011110111100000001101110';
elseif q==512
    CCSK_seq_str = '10100000110010000011000111111100011010100110101011010000101111100000011111010101101011000100010100000101011000001011011011011000101110100011001101000000000101111011101100100000011011111011111110001110101001100010000010111011000111010010110010000111100100101011110100011010001110000111100011110011001010100010101111001100011101101100100000110001000010000011010010001011100111111111010101000111000111100111100110011100100101111011100100100100111100101110011000110111101100111111011110110000000111110110010101001000';
end

n=log2(N);
related_variables_pth = fullfile(mainPath, 'related_variables');
pth3 = (fullfile(related_variables_pth, 'GF_arithm'));
fl_nm = ['arith_' num2str(q) '.mat'];
load(fullfile(pth3, fl_nm));
SNRs = 10.^(SNRs_db/10);
N0 = 1./(2*SNRs);
sigm = sqrt(N0);
words = (0:q-1);
p=log2(q);
eta0 = CCSK_seq_str' - 48;
etaq = zeros(q,q);
for i=0:q-1
    etaq(:,i+1)=circshift(eta0,-i);
end
etaqm = 1-2*etaq;

tic
for prf=1:length(sigm)
    cl_cn = cell(n,1);
    cl_vn = cell(n,1);
    ncl = zeros(n,1);
    nkrnl = zeros(n,1);
    for l1=n:-1:1
        a=nan(N/2,1);
        b=a;
        for t=0:N/2-1
            a(t+1)=2*t-mod(t,2^(l1-1))+1;
            b(t+1)=2^(l1-1)+2*t-mod(t,2^(l1-1))+1;
        end
        a1=reshape(a, 2^(l1-1),(N/2)/(2^(l1-1)));
        b1=reshape(b, 2^(l1-1),(N/2)/(2^(l1-1)));
        cl_cn{l1}=a1;
        cl_vn{l1}=b1;
        ncl(l1)=size(a1,2);
        nkrnl(l1) = N/ncl(l1)/2;
    end

    e_entr=cell(n,1);
    e_err_prob=cell(n,1);
    for l2=1:n
        e_entr{l2}=nan(n_obs,N);
        e_err_prob{l2}=nan(n_obs,N);
    end

    h1 = nan(N,n);

    for i = 1 : n_obs
        PP=nan(q,N, n+1);
        info_seq=randi([0 q-1], 1, N);
        [info_seq,x,m1, y]=gen_noisy_seq_ccsk(N,q,sigm(prf), add_mat, etaqm, info_seq);
        [P, HD_L]=cond_prob_ccsk(y,sigm(prf),etaq);
        PP(:,:,n+1)=P;
        for l2=n:-1:1
            for nc1=1:ncl(l2)
                for cl1=1:nkrnl(l2)
                    ii1 = cl_cn{l2}(cl1,nc1);
                    ii2 = cl_vn{l2}(cl1,nc1);
                    L1=squeeze(PP(:,ii1,l2+1));
                    L2=squeeze(PP(:,ii2,l2+1));
                    u0=m1(ii1, l2);
                    u1=m1(ii2, l2);
                    [Pr0, Pr1]= F1(L1, L2, add_mat,mul_mat, 1, u0);
                    PP(:,ii1,l2)=Pr0;
                    PP(:,ii2,l2)=Pr1;
                    e_entr{l2}(i, ii1)=cond_entropy(Pr0)/p;
                    e_entr{l2}(i, ii2)=cond_entropy(Pr1)/p;
                    e_err_prob{l2}(i, ii1) =1-Pr0(u0+1);
                    e_err_prob{l2}(i, ii2)=1-Pr1(u1+1);
                end
            end
        end
    end

    %%
    temp_emtrop = cell(n,1);
    temp_erp = cell(n,1);
    entrop = cell(n,1);
    erp = entrop;
    sorted_entropy = entrop;
    sorted_err_prob = entrop;
    ch_idx_entr = entrop;
    ch_idx_erp = entrop;
    for l=1:n
        temp_emtrop{l} = mean(e_entr{l},1);
        temp_erp{l} = mean(e_err_prob{l},1);
        for clst = 1:ncl(l)
            ii1=cl_cn{l}(:,clst);
            ii2=cl_vn{l}(:,clst);
            ent0 = mean(temp_emtrop{l}(ii1));
            ent1 = mean(temp_emtrop{l}(ii2));
            entrop{l}(2*clst-1:2*clst)=[ent0 ent1];
            erp0=mean(temp_erp{l}(ii1));
            erp1=mean(temp_erp{l}(ii2));
            erp{l}(2*clst-1:2*clst)=[erp0 erp1];
        end
    end

    for l=1:n
        [sorted_entropy{l}, ch_idx_entr{l}] = sort(entrop{l}, "ascend");
        [sorted_err_prob{l}, ch_idx_erp{l}] = sort(erp{l}, "ascend");
        if low2highreliab
            [sorted_entropy{l}, ch_idx_entr{l}] = sort((entrop{l}), "descend");
            [sorted_err_prob{l}, ch_idx_erp{l}] = sort((erp{l}), "descend");
        end

        ch_idx_entr{l} = ch_idx_entr{l}-1;
        ch_idx_erp{l}=ch_idx_erp{l}-1;
    end
    %%
    for l=1:n
        folderName = ['.\ccsk_reliab_seq\N' num2str(2^(n-l+1))];
        fileName = ['mat_N' num2str(2^(n-l+1)) '_GF' num2str(q) '_SNR' num2str(SNRs_db(prf), '%.3f') '.txt'];
        flnm = fullfile(folderName, fileName);


        if ~exist(folderName, 'dir')
            mkdir(folderName);
        end
        fileID = fopen(flnm, 'w');
        for i = 1:2^(n-l+1)
            fprintf(fileID, '%d ', ch_idx_entr{l}(i));
        end
        fprintf(fileID, '\n\n');
        for i = 1:2^(n-l+1)
            fprintf(fileID, '%d ', ch_idx_erp{l}(i));
        end

        fprintf(fileID, '\n\n');
        for i = 1:2^(n-l+1)
            fprintf(fileID, '%.12f ', entrop{l}(i));
        end
        fprintf(fileID, '\n\n');
        for i = 1:2^(n-l+1)
            fprintf(fileID, '%.12f ', erp{l}(i));
        end
        fprintf(fileID, '\n\nN_obs: %d',n_obs*(2^(l-1)) );
        fprintf(fileID, '\n\nFile Format:\n');
        if low2highreliab
            fprintf(fileID, '1- least to most reliable channels (taken entropies)\n');
            fprintf(fileID, '2- least to most reliable channels (taken error probabilities)\n');
        else
            fprintf(fileID, '1- most to least reliable channels (taken entropies)\n');
            fprintf(fileID, '2- most to least reliable channels (taken error probabilities)\n');
        end
        fprintf(fileID, '3- channels entropies\n');
        fprintf(fileID, '4- channels error probabilities\n');
        fprintf(fileID, '5- number of MontCarlo simulations');

        fclose(fileID);
    end

end
toc


