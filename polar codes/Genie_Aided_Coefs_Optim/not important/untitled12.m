clear
close all
% rng(1);
hamming_distance = @(a, b) arrayfun(@(x) sum(dec2bin(bitxor(a, x)) == '1'), b);

pth1 = (fullfile(pwd, 'related_functions'));
addpath(pth1);
mainPath = pwd;
q=64;
N=2;
opt_entr=1;
n_obs=5000;

n=log2(N);
p=log2(q);
SNRs_db = 0;
SNRs = 10.^(SNRs_db/10);
N0 = 1./(SNRs);
sigm = sqrt(N0);
words = (0:q-1);
alph_bin =  fliplr(dec2bin(words, p) - 48);
alph_bin_mod = (-1).^alph_bin;
h01=1;

related_variables_pth = fullfile(mainPath, 'related_variables');
pth3 = (fullfile(related_variables_pth, 'GF_arithm'));
fl_nm = ['arith_' num2str(q) '.mat'];
load(fullfile(pth3, fl_nm));


t1d=[0	1	2	4	8	16	32	3	5	9	17	33	6	10	18	34	12	20	36	24	40	48	7	11	19	35	13	21	37	25	41	49	14	22	38	26	42	50	28	44	52	56	15	23	39	27	43	51	29	45	53	57	30	46	54	58	60	31	47	55	59	61	62	63];
[word1,ii1]=sort(t1d);
t2d=[0	31	47	55	59	61	62	48	40	36	34	33	24	20	18	17	12	10	9	6	5	3	7	11	13	14	19	21	22	25	26	28	35	37	38	41	42	44	49	50	52	56	60	58	57	54	53	51	46	45	43	39	30	29	27	23	15	1	2	4	8	16	32	63];
% t2d=[0	31	47	48	55	40	24	7	59	36	20	11	12	19	35	60	61	34	18	13	10	21	37	58	6	25	41	54	49	46	30	1	62	33	17	14	9	22	38	57	5	26	42	53	50	45	29	2	3	28	44	51	52	43	27	4	56	39	23	8	15	16	32	63];
word1_p=t2d(ii1);
% u = randi([0 q-1],1,2);
u=[10 18];
X = [add_mat(u(1)+1, u(2)+1) word1_p(u(2)+1)];

word1_bin = dec2bin(word1,p)-48;

en0=zeros(n_obs,1);
en1=zeros(n_obs,1);

Delta_max=0;

P1=zeros(q,n_obs);
P2=zeros(q,n_obs);

Pr0m=P1;
Pr1m=P1;
for no=1:n_obs
    x=X(1);
    etax = alph_bin_mod(x+1, :);
    nse = sigm*randn(size(etax));
    y = etax + nse;
    [P1(:,no), ~]=cond_prob(y, sigm,alph_bin);
    en0(no) =cond_entropy(P1(:,no));

    x=X(2);
    etax = alph_bin_mod(x+1, :);
    nse = sigm*randn(size(etax));
    y = etax + nse;
    [P2(:,no), ~]=cond_prob(y, sigm,alph_bin);
    en1(no) =cond_entropy(P2(:,no));
    [Pr0m(:,no), Pr1m(:,no)]= F1_permutation(P1(:,no), P2(:,no), add_mat,word1_p, u(1));

end

%%
Pm1=mean(P1,2);
Pm2=mean(P2,2);

% [Pms1, idx1]=sort(Pm1,'descend');
% idx1=idx1-1;
% % idx1_b=de2bi(idx1,p);
% lst1 = idx1(2:nchoosek(6,1)+nchoosek(6,2)+1)';
% d = hamming_distance(X(1), idx1);

% [Pms2, idx2]=sort(Pm2,'descend');
% idx2=idx2-1;
% % idx_b1=de2bi(idx2,p);
% lst2 = idx2(2:nchoosek(6,1)+nchoosek(6,2)+1)';
% 
% u0=add_mat(X(1)+1, X(2)+1);
% lst3=add_mat(u0+1, lst1+1);


% 
% intsct1=intersect(lst3, lst2);
% 
% Pm3=zeros(q,1);
% Pm5=Pm3;
% Pm3(u0+1)=1;
% for a=0:q-1
%     for b = 0:q-1
%         u=add_mat(1+a, 1+b);
%         Pm5(u+1)= Pm5(u+1)+Pm1(a+1)*Pm3(b+1);
%     end
% end




% Pm4=Pm5.*Pm2;
% Pm5=Pm5/sum(Pm5);
% Pm5=Pm5/sum(Pm5);
% Pm4=Pm4/sum(Pm4);
% nn=63;
% [Pm4_s, Pm4_s_i]=sort(Pm4);
% il1 = Pm4_s_i(1:nn)-1;
% ih1 = flip(Pm4_s_i(q-nn:q-1)-1);
% word1_p1=(0:q-1)';

% word1_p = mul_mat(word1+1, 50+1);


[Pr0, Pr1]= F1_permutation(Pm1, Pm2, add_mat,word1_p, u(1));
% plot(Pm2)
hold on
plot(Pr1)
hold off
%%
% u(1)=add_mat(X(1)+1,X(2)+1);
% u(2)=X(2);
for no=1:n_obs
    u=randi([0 q-1],1,2);
    x=u;
    x(1) = add_mat(1+u(1), 1+u(2));
    x(2)=word1_p(u(2)+1);
    etax = alph_bin_mod(x+1, :);
    nse = sigm*randn(size(etax));
    y = etax + nse;
    [P, HD_L]=cond_prob(y, sigm,alph_bin);
    L1=P(:,1);
    L2=P(:,2);
    u0=u(1);
    [Pr0, Pr1]= F1_permutation(L1, L2, add_mat,word1_p, u0);
    [~,aa]=max((Pr0));
    en0(no) =cond_entropy(Pr0);
    en1(no)=cond_entropy(Pr1);
end

mne0=p-mean(en0);
mne1=p-mean(en1);
Delta=abs(mne0-mne1);
if(Delta>Delta_max)
    Delta_max=Delta;
end
Delta_max
