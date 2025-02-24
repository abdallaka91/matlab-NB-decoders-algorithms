clear
% rng(1);
pth1 = (fullfile(pwd, 'related_functions'));
addpath(pth1);
mainPath = pwd;
q=16;
N=2;
opt_entr=1;
n_obs=1000;


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

word1= 0:q-1;
word1_bin = dec2bin(word1,p)-48;

en0=zeros(n_obs,1);
en1=zeros(n_obs,1);

h1=63;
mne0=zeros(h1,1);
mne1=mne0;
Delta=mne0;
Delta_max=0;

u = [13 5];
t1d=[0 1 9 13 15 14 6 2 3 8 5 11 12 7 10 4];
[t1d,ii1]=sort(t1d);
t2d=[0 7 9 2 15 8 6 13 10 14 12 4 5 1 3 11];
word1_p=t2d(ii1);

for no=1:n_obs
    %         u=randi([0 q-1],1,N);
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
    %         [Pr0, Pr1]= F1(L1, L2, add_mat,mul_mat, h1, u0);
    [~,aa]=max((Pr0));
    en0(no) =cond_entropy(Pr0);
    en1(no)=cond_entropy(Pr1);

end

mne0(h1)=mean(en0);
mne1(h1)=mean(en1);
Delta(h1)=abs(mne0(h1)-mne1(h1));
if(Delta(h1)>Delta_max)
    Delta_max=Delta(h1);
    h11=h1;
end

