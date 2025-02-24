clear
% rng(1);

pth1 = (fullfile(pwd, 'related_functions'));
addpath(pth1);
mainPath = pwd;
opt_entr=1;
n_obs=1000;
q=64;
p=log2(q);
N=2;
n=log2(N);
related_variables_pth = fullfile(mainPath, 'related_variables');
pth3 = (fullfile(related_variables_pth, 'GF_arithm'));
fl_nm = ['arith_' num2str(q) '.mat'];
load(fullfile(pth3, fl_nm));

SNRs_db = 0;
SNRs = 10.^(SNRs_db/10);
N0 = 1./(SNRs);
sigma = sqrt(N0);

word1= 0:q-1;
word1_bin = dec2bin(word1,p)-48;

word1_p=[54	49	7	26	42	4	61	9	0	29	43	44	19	39	30	48	28	40	17	63	15	18	36	35	37	11	50	6	57	62	8	21	45	3	32	20	24	31	51	46	14	58	25	55	52	41	5	2	59	38	10	13	1	53	22	56	23	16	60	33	34	12	47	27];
% 56	1	20	55	31	43	2	44	51	29	46	26	37	6	9	48	45	14	27	34	54	24	49	5	0	52	7	41	42	19	28	63	10	36	13	57	17	50	39	30	22	47	32	3	60	8	59	21	35	23	62	16	4	61	40	11	25	58	53	12	15	33	18	38
t1d=word1;
t2d=word1_p;

alph_bin =  fliplr(dec2bin(word1, p) - 48);
alph_bin_mod = (-1).^alph_bin;

[sp_frwd, sp_bck, hmdst1, hmdst2]=hamm_sp_perm(t1d, t2d);

sigm = sigma(1);

mne0=cell(n,1);
mne1=mne0;
e0=mne0;
e1=mne1;
Deltas = mne1;

mx1=0;

ttt=zeros(1,q);
tttt=ttt;
for n0 = 1:n

    mne0{n-n0+1}=nan(1,2^(n0-1));
    mne1{n-n0+1}=nan(1,2^(n0-1));
    e0{n-n0+1}=nan(1,n_obs, 2^(n0-1));
    e1{n-n0+1}=nan(1,n_obs, 2^(n0-1));

    N1=2^n0;
    cl_cn = cell(n0,1);
    cl_vn = cell(n0,1);
    ncl = zeros(n0,1);
    for l1=n0:-1:1
        a=nan(N1/2,1);
        b=a;
        for t=0:N1/2-1
            a(t+1)=2*t-mod(t,2^(l1-1))+1;
            b(t+1)=2^(l1-1)+2*t-mod(t,2^(l1-1))+1;
        end
        a1=reshape(a, 2^(l1-1),(N1/2)/(2^(l1-1)));
        b1=reshape(b, 2^(l1-1),(N1/2)/(2^(l1-1)));
        cl_cn{l1}=a1;
        cl_vn{l1}=b1;
        ncl(l1)=size(a1,2);
    end
    PP=nan(q,N1, n0+1);
    for i = 1 : n_obs

%         info_seq=randi([0 q-1], 1, N1);
        info_seq=[10 18];
        i2=cl_vn{1}(:);
        [info_seq,x,m1, y]=encode_and_nois_permut(N1,q,sigm,word1_p, add_mat, alph_bin_mod, 1, info_seq);
        for l2=n0:-1:1
            for nc1=1:ncl(l2)
                for cl1=1:size(cl_cn{l2},1)
                    ii1 = cl_cn{l2}(cl1, nc1);
                    ii2 = cl_vn{l2}(cl1, nc1);
                    if l2==n0
                        [P, HD_L]=cond_prob(y([ii1 ii2],:),sigm,alph_bin);
                        PP(:,[ii1 ii2],n0+1)=P;
                    end
                    L1=squeeze(PP(:,ii1,l2+1));
                    L2=squeeze(PP(:,ii2,l2+1));
                    u0=m1(ii1, l2);
                    u1=m1(ii2, l2);
                    if l2>l1
                        [Pr0, Pr1]= F1_permutation(L1, L2, add_mat,word1_p, u0);
                        PP(:,ii1,l2)=Pr0;
                        PP(:,ii2,l2)=Pr1;
                    else
                        L1=squeeze(PP(:,ii1,l1+1));
                        L2=squeeze(PP(:,ii2,l1+1));
                        [Pr0, Pr1]= F1_permutation(L1, L2, add_mat,word1_p, u0);
                        ttt=ttt+Pr1;
                            tttt=tttt+L2;
                        if ~opt_entr
                            e0{n-n0+1}(i, nc1) =1-Pr0(u0+1);
                            e1{n-n0+1}(i, nc1)=1-Pr1(u1+1);
                        else
                            e0{n-n0+1}(i, nc1) =cond_entropy(Pr0);
                            e1{n-n0+1}(i, nc1)=cond_entropy(Pr1);
                        end
                    end
                end
            end
        end
    end

    for nc1=1:2^(n0-1)
        mne0{n-n0+1}( nc1)=log2(q)-mean(squeeze(e0{n-n0+1}(:, nc1)));
        mne1{n-n0+1}( nc1)=log2(q)-mean(squeeze(e1{n-n0+1}(:, nc1)));
        Deltas{n-n0+1}( nc1) = abs(mne0{n-n0+1}( nc1)- mne1{n-n0+1}( nc1));
    end

end
%     jj
%     mne0
%     mne1
% if Deltas{1}>mx1
%     mx1=Deltas{1};
%     word1_pb=word1_p;
% end
%     deltas1(jj)=Deltas{1};
% end
% Deltas{1}
