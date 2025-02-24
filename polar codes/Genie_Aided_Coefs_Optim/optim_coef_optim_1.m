clear
% rng(1);
pth1 = (fullfile(pwd, 'related_functions'));
addpath(pth1);
mainPath = pwd;
q=64;
p=log2(q);
N=256;
n=log2(N);
related_variables_pth = fullfile(mainPath, 'related_variables');
pth3 = (fullfile(related_variables_pth, 'GF_arithm'));
fl_nm = ['arith_' num2str(q) '.mat'];
load(fullfile(pth3, fl_nm));

SNRs_db = 0;
SNRs = 10.^(SNRs_db/10);
N0 = 1./(SNRs);
sigma = sqrt(N0);
words = (0:q-1);
alph_bin =  fliplr(dec2bin(words, p) - 48);
alph_bin_mod = (-1).^alph_bin;
h01=1;

sigm = sigma(1);
opt_entr=1;
n_obs=1000;

coefs=cell(n,1);
for i=1:n
    coefs{i}=64+2^i+(0:2^(i-1)-1);
end

mne0=cell(n,1);
mne1=mne0;
e0=mne0;
e1=mne1;
Deltas=mne0;
tic

coefsp={[                                                                                     25]
    [                                                                                   3 18]
    [                                                                            44 42 13 23]
    [                                                                35 59 21 41 39 61 34 34]
    [                                         42 19 42 36 59 8 61 41 61 41 29 50 18 54 33 39]
    [32 25 59 36 44 43 3 7 5 1 17 7 52 7 48 46 18 29 50 24 57 46 15 2 46 46 56 12 56 34 51 1]
    [58 14 45 9 16 33 38 62 2 4 16 26 27 60 34 60 40 14 3 52 5 55 4 29 35 62 60 42 4 31 27 18 29 14 45 11 30 46 17 27 43 42 59 16 23 43 15 1 14 23 55 48 24 28 31 1 9 4 38 1 42 1 1 1]
    };

nnn=length(coefsp);

coefs(1:nnn)=coefsp(1:nnn);

for n0 = nnn+1:n

    mne0{n-n0+1}=nan(q-1,2^(n0-1));
    mne1{n-n0+1}=nan(q-1,2^(n0-1));
    e0{n-n0+1}=nan(q-1,n_obs, 2^(n0-1));
    e1{n-n0+1}=nan(q-1,n_obs, 2^(n0-1));

    N1=2^n0;
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
    for i=1:n0-1
        for j=1:2^(i-1)
            h1(cl_vn{n0-i+1}(:,j),n0-i+1)=coefs{i}(j);
        end
    end
    for i = 1 : n_obs
        PP=nan(q,N1, n0+1);
        info_seq=randi([0 q-1], 1, N1);
        for h0 = 1:q-1
            i2=cl_vn{1}(:);
            h1(i2,1)=h0;
            [info_seq,x,m1, y]=gen_noisy_seq(N1,q,sigm,h1, add_mat, mul_mat, alph_bin_mod, 1, info_seq);
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
                        if l2>1
                            h2 = h1(ii2,l2);
                            [Pr0, Pr1]= F1(L1, L2, add_mat,mul_mat, h2, u0);
                            PP(:,ii1,l2)=Pr0;
                            PP(:,ii2,l2)=Pr1;
                        else
                            L1=squeeze(PP(:,ii1,2));
                            L2=squeeze(PP(:,ii2,2));
                            [Pr0, Pr1]= F1(L1, L2, add_mat,mul_mat, h0, u0);
                            if ~opt_entr
                                e0{n-n0+1}(h0,i, nc1) =1-Pr0(u0+1);
                                e1{n-n0+1}(h0,i, nc1)=1-Pr1(u1+1);
                            else
                                e0{n-n0+1}(h0,i, nc1) =p-cond_entropy(Pr0);
                                e1{n-n0+1}(h0,i, nc1)=p-cond_entropy(Pr1);
                            end
                        end
                    end
                end
            end
        end
    end
    for h0=1:q-1
        for nc1=1:2^(n0-1)
            mne0{n-n0+1}(h0, nc1)=mean(squeeze(e0{n-n0+1}(h0,:, nc1)));
            mne1{n-n0+1}(h0, nc1)=mean(squeeze(e1{n-n0+1}(h0,:, nc1)));
            Deltas{n-n0+1}(h0, nc1)=abs(mne0{n-n0+1}(h0, nc1)-mne1{n-n0+1}(h0, nc1));
        end
    end
    for nc1=1:2^(n0-1)
        aa=abs(mne0{n-n0+1}(:,nc1)-mne1{n-n0+1}(:,nc1));
        [~,ii]=sort(aa,'descend');
        coefs{n0}(nc1)=ii(1);
    end

end
toc
%%
for l1=n:-1:1
    b=nan(N/2,1);
    for t=0:N/2-1
        b(t+1)=2^(l1-1)+2*t-mod(t,2^(l1-1))+1;
    end
    b1=reshape(b, 2^(l1-1),(N/2)/(2^(l1-1)));
    cl_vn{l1}=b1;
    for i =1 : size(cl_vn{l1},2)
        s1=cl_vn{l1}(:,i);
        h1(s1,l1)=coefs{n-l1+1}(i);
    end
end

entropy_prob = zeros(1,N);
if opt_entr==1
    opt1='Entropy';
    for k=1:N/2
        ent0=mne0{1,1}(h1(2*k,1),k);
        ent1=mne1{1,1}(h1(2*k,1),k);
        entropy_prob(2*k-1:2*k)=[ent0 ent1];
    end
    [sorted_entropy_prob, ch_idx] = sort(entropy_prob, "descend");
else
    opt1='error_prob';
    for k=1:N/2
        ent0=mne0{1,1}(h1(2*k,1),k);
        ent1=mne1{1,1}(h1(2*k,1),k);
        entropy_prob(2*k-1:2*k)=[ent0 ent1];
    end
    [sorted_entropy_prob, ch_idx] = sort(entropy_prob, "ascend");
end
Deltas1=abs(entropy_prob(1:2:end)-entropy_prob(2:2:end));
for l1=n:-1:1
    b=nan(N/2,1);
    for t=0:N/2-1
        b(t+1)=2^(l1-1)+2*t-mod(t,2^(l1-1))+1;
    end
    b1=reshape(b, 2^(l1-1),(N/2)/(2^(l1-1)));
    cl_vn{l1}=b1;
    for i =1 : size(cl_vn{l1},2)
        s1=cl_vn{l1}(:,i);
        h1(s1,l1)=coefs{n-l1+1}(i);
    end
end

save([pwd '\optimized_coefs\N' num2str(N), '_q' num2str(q) '_snr' num2str(SNRs_db+0.01,3) 'dB_observ' num2str(n_obs) '_' opt1 '.mat'],'h1','sorted_entropy_prob','ch_idx', 'opt1', 'q' )
save([pwd '\optimized_coefs\N' num2str(N), '_q' num2str(q) '_snr' num2str(SNRs_db+0.01,3) 'dB_observ' num2str(n_obs) '_' opt1 '_all_param.mat'])
