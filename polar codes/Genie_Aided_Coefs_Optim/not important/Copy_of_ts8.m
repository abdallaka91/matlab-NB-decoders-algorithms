clear
% rng(1);
pth1 = (fullfile(pwd, 'related_functions'));
addpath(pth1);
mainPath = pwd;

SNRs_db = 0;
n_obs=2000;
q=64;
p=log2(q);
N=4;
opt_entr=1;

% data1=load(['valen_0dB_N' num2str(8) '.mat']);
data1=load('C:\Users\Abdallah Abdallah\Documents\Personals\matlab-NB-decoders-algorithms\polar codes\Genie_Aided_Coefs_Optim\optimized_coefs\N4_q64_snr0.01dB_observ1000_Entropy.mat');
h1=data1.h1;

n=log2(N);
related_variables_pth = fullfile(mainPath, 'related_variables');
pth3 = (fullfile(related_variables_pth, 'GF_arithm'));
fl_nm = ['arith_' num2str(q) '.mat'];
load(fullfile(pth3, fl_nm));


SNRs = 10.^(SNRs_db/10);
N0 = 1./(SNRs);
sigma = sqrt(N0);
words = (0:q-1);
alph_bin =  fliplr(dec2bin(words, p) - 48);
alph_bin_mod = (-1).^alph_bin;
h01=1;

sigm = sigma(1);



mne0=cell(n,1);
mne1=mne0;
e0=mne0;
e1=mne1;
Deltas=mne0;
n0=n;
N1=2^n0;

iii=0;
% valentin coefs at 3 dB [38 ; 35 35;  41 53 29 41]
% coefs1={randi([1 q-1],1,1), randi([1 q-1],1,2), randi([1 q-1],1,4)};
% coefs1={38 ; [35 35];  [41 53 29 41]};
coefs1=cell(n,1);

for l1=n0:-1:1
    b=nan(N1/2,1);
    for t=0:N1/2-1
        b(t+1)=2^(l1-1)+2*t-mod(t,2^(l1-1))+1;
    end
    b1=reshape(b, 2^(l1-1),(N1/2)/(2^(l1-1)));
    cl_vn{l1}=b1;
    coefs1{n0-l1+1}=(h1(cl_vn{l1}(1,:),l1))';
end

coefs2=coefs1;
cnd1=false;
coefs_lst = nan(100,sum(cellfun(@numel, coefs1)));
chn_MI_lst=nan(100,N);
i9=1;

cnd2=false;
mx=0;
while iii<10
    if mod(iii,1)==0 && iii>0
        if cnd2==false
            for i10=1:n
                coefs2{i10}=randi([1 q-1],1, 2^(i10-1));
            end
        else
            cnd2=false;
        end
    end
    

    coefs=coefs2;
    for i4 = 1:numel(coefs1)
        for j4 = 1:numel(coefs1{i4})
            
            mx1=0;
            for v4=0:q-2
                values = coefs1{i4}; 
                cc1=values(j4);
                values(j4)=mod((cc1+v4)-1,q-1)+1;
                coefs{i4}=values;
                h1 = nan(N1,n0);
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
                for i=1:n0
                    for j=1:2^(i-1)
                        h1(cl_vn{n0-i+1}(:,j),n0-i+1)=coefs{i}(j);
                    end
                end

                for i=1:n0
                    e0{i}=nan(n_obs*2^(i-1), ncl(i));
                    e1{i}=e0{i};
                    mne0{i}=nan(1, ncl(i));
                    mne1{i}=nan(1, ncl(i));
                end
                for i = 1 : n_obs
                    PP=nan(q,N1, n0+1);
                    info_seq=randi([0 q-1], 1, N1);
                    [info_seq,x,m1, y]=gen_noisy_seq(N1,q,sigm,h1, add_mat, mul_mat, alph_bin_mod, 1, info_seq);
                    for l2=n0:-1:1
                        for nc1=1:ncl(l2)
                            temp0 = nan(1,size(cl_cn{l2},1));
                            temp1 = nan(1,size(cl_cn{l2},1));
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
                                h2 = h1(ii2,l2);
                                [Pr0, Pr1]= F1(L1, L2, add_mat,mul_mat, h2, u0);
                                PP(:,ii1,l2)=Pr0;
                                PP(:,ii2,l2)=Pr1;

                                if ~opt_entr
                                    temp0(cl1) =1-Pr0(u0+1);
                                    temp1(cl1)=1-Pr1(u1+1);
                                else
                                    temp0(cl1) =cond_entropy(Pr0);
                                    temp1(cl1)=cond_entropy(Pr1);
                                end

                            end
                            e0{l2}((i-1)*length(temp0)+1:i*length(temp0),nc1)=temp0;
                            e1{l2}((i-1)*length(temp0)+1:i*length(temp0), nc1)=temp1;
                        end
                    end
                end
                %
                for l2=n0:-1:1
                    for nc1=1:ncl(l2)
                        mne0{l2}(1,nc1)=log2(q)-mean(e0{l2}(:,nc1));
                        mne1{l2}(1,nc1)=log2(q)-mean(e1{l2}(:,nc1));
                        Deltas{l2}(1,nc1) = abs(mne0{l2}(1,nc1)-mne1{l2}(1,nc1));
                    end
                end
                %      coefs
                %     Deltas
                chn_MI = [mne0{1};mne1{1}];
                chn_MI=chn_MI(:);
                chn_MI=chn_MI';
                chn_MI_srt=sort(chn_MI);
                mn=mean(chn_MI_srt);
                sm1=sum(abs((chn_MI_srt-mn)));
                if sm1>mx1
                    mx1=sm1;
                    coefs2{i4}(j4)=values(j4);
                end

                if sm1>mx
                    cnd2=true;
                    mn1=mn;
                    mx=sm1;
                    coefs1=coefs;
                    chn_MI1=chn_MI;
                    chn_MI_lst(i9,:)=chn_MI1;
                    i9=i9+1;
                    1/[1 2]
                end
            end
            coefs{i4}(j4)=coefs2{i4}(j4);
        end
    end
    iii=iii+1;
%     if iii==1
%         coefs11=coefs1;
%         cnd1=false;
%     elseif isequal(coefs11, coefs1)
%         cnd1=true;
%     else
%         coefs11=coefs1;
%     end
end