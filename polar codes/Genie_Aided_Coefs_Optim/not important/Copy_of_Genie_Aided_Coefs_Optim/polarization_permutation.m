clear
% rng(1);
pth1 = (fullfile(pwd, 'related_functions'));
addpath(pth1);
mainPath = pwd;
q=64;
N=2;
opt_entr=0;
n_obs=10000;


n=log2(N);
p=log2(q);
SNRs_db = 1;
SNRs = 10.^(SNRs_db/10);
N0 = 1./(SNRs);
sigm = sqrt(N0);
words = (0:q-1);
alph_bin =  fliplr(dec2bin(words, p) - 48);
alph_bin_mod = (-1).^alph_bin;


related_variables_pth = fullfile(mainPath, 'related_variables');
pth3 = (fullfile(related_variables_pth, 'GF_arithm'));
fl_nm = ['arith_' num2str(q) '.mat'];
load(fullfile(pth3, fl_nm));

word1= 0:q-1;
word1_bin = dec2bin(word1,p)-48;

h1=63;
mne0=zeros(h1,1);
mne1=mne0;
Delta=mne0;
Delta_max=0;

U = [0 7];

per0=0:q-1;
% per1=per0;
% per0=mul_mat(50+1,:);
% per0=[0	1	2	15	4	23	27	7	8	29	30	11	39	13	14	3	16	43	45	19	46	21	22	5	51	25	26	6	28	9	10	31	32	53	54	35	57	37	38	12	58	41	42	17	44	18	20	47	60	49	50	24	52	33	34	55	56	36	40	59	48	61	62	63];
% per1=[0	8	16	24	25	40	48	7	1	9	10	26	12	28	32	57	2	17	18	41	20	42	44	58	3	4	11	27	13	29	30	49	14	33	34	50	36	52	59	60	5	19	21	43	22	45	46	61	6	31	35	51	37	53	54	62	56	15	23	38	39	47	55	63];%optim for 7
per1=[56	1	20	55	31	43	2	44	51	29	46	26	37	6	9	48	45	14	27	34	54	24	49	5	0	52	7	41	42	19	28	63	10	36	13	57	17	50	39	30	22	47	32	3	60	8	59	21	35	23	62	16	4	61	40	11	25	58	53	12	15	33	18	38];
% per1=mul_mat(50+1,:);
per0_inv(per0+1)= 0:length(per0)-1;
per1_inv(per1+1)= 0:length(per1)-1;

P0=zeros(q,n_obs);
P1=zeros(q,n_obs);

Pr0m=P0;
Pr1m=P1;

en0=zeros(n_obs,1);
en1=zeros(n_obs,1);

P=nan(q,2,n_obs);

for no=1:n_obs

    x0=U;
    x0(1) = add_mat(1+U(1), 1+U(2));
    x0(2)=U(2);

    x(1)=per0(x0(1)+1);
    x(2)=per1(x0(2)+1);
    etax = alph_bin_mod(x+1, :);
    nse = sigm*randn(size(etax));
    y = etax + nse;
    [P(:,:,no), HD_L]=cond_prob(y, sigm,alph_bin);
    L0=squeeze(P(:,1,no));
    L1=squeeze(P(:,2,no));
    u0=U(1);
    [Pr0m(:,no), Pr1m(:,no)]= F2_permutation(L0, L1, add_mat,per0,per1, u0);
        en0(no) =cond_entropy(Pr0m(:,no));
        en1(no) =cond_entropy(Pr1m(:,no));
%     en0(no) =1-(Pr0m(U(1),no));
%     en1(no) =1-(Pr1m(U(2),no));

end

Pm0=mean(Pr0m,2);
Pm1=mean(Pr1m,2);

L0a=mean(P(:,1,:),3);
L1a=mean(P(:,2,:),3);

[Pr0a, Pr1a]= F2_permutation(L0a, L1a, add_mat,per0,per1, u0);

L0a_p(per0+1) = L0a;
L1a_p(per1+1) = L1a;

L2a=zeros(1,q);
L4a=L2a;
L2a(u0+1)=1;

for a=0:q-1
    for b = 0:q-1
        u=add_mat(1+a, 1+b);
        L4a(u+1)= L4a(u+1)+L0a_p(a+1)*L2a(b+1);
    end
end

mean(en0)-mean(en1)

% figure;
% stem(words,L1a_p/sum(L1a_p),'k+')
% hold on
% stem(words,L4a/sum(L4a),'r.')

% figure;
% stem(words,Pr1a/sum(Pr1a),'m+')
% hold on
% stem(words,L1a_p.*L4a/sum(L1a_p.*L4a),'o')
