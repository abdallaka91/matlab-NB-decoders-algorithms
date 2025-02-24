clear
% rng(1);
hamming_distance = @(a, b) arrayfun(@(x) sum(dec2bin(bitxor(a, x)) == '1'), b);

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
u = [13 5];
t1d=[0	1	2	4	8	16	32	3	5	9	17	33	6	10	18	34	12	20	36	24	40	48	7	11	19	35	13	21	37	25	41	49	14	22	38	26	42	50	28	44	52	56	15	23	39	27	43	51	29	45	53	57	30	46	54	58	60	31	47	55	59	61	62	63];
[word1,ii1]=sort(t1d);
t2d=[0	31	47	55	59	61	62	48	40	36	34	33	24	20	18	17	12	10	9	6	5	3	7	11	13	14	19	21	22	25	26	28	35	37	38	41	42	44	49	50	52	56	60	58	57	54	53	51	46	45	43	39	30	29	27	23	15	1	2	4	8	16	32	63];
word1_p=t2d(ii1);
% word1_p=[0	42	44	6	47	5	3	41	15	37	35	9	32	10	12	38	31	53	51	25	48	26	28	54	16	58	60	22	63	21	19	57	23	61	59	17	56	18	20	62	24	50	52	30	55	29	27	49	8	34	36	14	39	13	11	33	7	45	43	1	40	2	4	46];
alph_bin =  fliplr(dec2bin(word1, p) - 48);
alph_bin_mod = (-1).^alph_bin;

cmbs = nchoosek(1:q, 2);
lc=size(cmbs,1);

hmdst1=zeros(lc,1);
for i=1:lc
    hmdst1(i)=hamming_distance(t1d(cmbs(i,1)), t1d(cmbs(i,2)));
end

hmdst2=zeros(lc,1);
for i=1:lc
    hmdst2(i)=hamming_distance(t2d(cmbs(i,1)), t2d(cmbs(i,2)));
end

dst1=1;

vv=[hmdst1 hmdst2];
idx1=hmdst1==1;
s1=sum(hmdst2(idx1));
idx2=hmdst1==2;
s2=sum(hmdst2(idx2));
vv1=[hmdst1(idx1) hmdst2(idx1)];
vv2=[hmdst1(idx2) hmdst2(idx2)];

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
Deltas{1}
