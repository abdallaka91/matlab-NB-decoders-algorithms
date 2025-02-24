clear
% rng(1);
pth1 = (fullfile(pwd, 'related_functions'));
addpath(pth1);
mainPath = pwd;
q=64;
p=log2(q);
N=64;
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

h1 = nan(N,n);
sigm = sigma(1);

%%
opt_entr=1;
obsrv_cnt=500;
PP=nan(q,N, n+1);
cl_cn = cell(n,1);
cl_vn = cell(n,1);
ncl = zeros(n,1);

for l1=n:-1:1
    a=nan(N/2,1);
    b=a;
    for t=0:N/2-1
        a(t+1)=2*t-mod(t,2^(l1-1))+1;
        b(t+1)=2^(l1-1)+2*t-mod(t,2^(l1-1))+1;
    end
    a1=reshape(a, 2^(l1-1),(N/2)/(2^(l1-1)));
    b1=reshape(b, 2^(l1-1),(N/2)/(2^(l1-1)));
    cl_cn{l1,1}=a1;
    cl_vn{l1,1}=b1;
    h1(b1,l1)=2;
    ncl(l1)=size(a1,2);
end

cpl=cell(n,n);
for l1=n:-1:1
    for l2=n:-1:l1
        cn1=downsample(cl_cn{l2}, size(cl_cn{l2},1)/(2^(l2-l1)));
        vn1=downsample(cl_vn{l2}, size(cl_vn{l2},1)/(2^(l2-l1)));
        cn1=cn1(:);
        vn1=vn1(:);
        cpl{l1,l2}=[cn1 vn1];
    end
end

mne0=cell(n,1);
mne1=mne0;
e0=mne0;
e1=mne1;

tic
for l1=n:-1:1
    mne0{l1,1}=nan(q-1,ncl(l1));
    mne1{l1,1}=nan(q-1,ncl(l1));
    e0{l1,1}=nan(q-1,obsrv_cnt, ncl(l1));
    e1{l1,1}=nan(q-1,obsrv_cnt, ncl(l1));
    %         id1=cl_vn{l1}(:,nc);
    id1 = find(~isnan(h1(:,l1)));
    for h0=1:q-1
        i=1;
        h1(id1,l1)=h0;
        info_seq=randi([0 q-1], 1, N);
        for i=1:obsrv_cnt
            PP=nan(q,N, n+1);
            [info_seq,x,m1, y]=gen_noisy_seq(N,q,sigm,h1, add_mat, mul_mat, alph_bin_mod, l1, info_seq);

            for l2=n:-1:l1
                ipx=cpl{l1,l2};
                for nc1=1:ncl(l1)
                    ii1 = ipx(nc1,1);
                    ii2 = ipx(nc1,2);
                    if l2==n
                        [P, HD_L]=cond_prob(y([ii1 ii2],:),sigm,alph_bin);
                        PP(:,[ii1 ii2],n+1)=P;
                    end
                    L1=squeeze(PP(:,ii1,l2+1));
                    L2=squeeze(PP(:,ii2,l2+1));
                    u0=m1(ii1, l2);
                    u1=m1(ii2, l2);
                    if l2>l1
                        h2 = h1(ii2,l2);
                    else
                        h2=h0;
                    end
                    [Pr0, Pr1]= F1(L1, L2, add_mat,mul_mat, h2, u0);
                    PP(:,ii1,l2)=Pr0;
                    PP(:,ii2,l2)=Pr1;
                end
            end
            for nc1=1:ncl(l1)
                ii1 = ipx(nc1,1);
                ii2=ipx(nc1,2);
                u0=m1(ii1, l2);
                u1=m1(ii2, l2);
                Pr0=PP(:,ii1,l1);
                Pr1=PP(:,ii2,l1);
                if ~opt_entr
                    e0{l1,1}(h0,i, nc1) =1-Pr0(u0+1);
                    e1{l1,1}(h0,i, nc1)=1-Pr1(u1+1);
                else
                    e0{l1,1}(h0,i, nc1) =cond_entropy(Pr0);
                    e1{l1,1}(h0,i, nc1)=cond_entropy(Pr1);
                end
                mne0{l1,1}(h0, nc1)=mean(squeeze(e0{l1,1}(h0,:, nc1)));
                mne1{l1,1}(h0, nc1)=mean(squeeze(e1{l1,1}(h0,:, nc1)));
            end
        end
        for nc1=1:ncl(l1)
            mne0{l1,1}(h0, nc1)=mean(squeeze(e0{l1,1}(h0,:, nc1)));
            mne1{l1,1}(h0, nc1)=mean(squeeze(e1{l1,1}(h0,:, nc1)));
        end
    end
    for nc1=1:ncl(l1)

        aa=abs(mne0{l1,1}(:,nc1)-mne1{l1,1}(:,nc1));
        [~,ii]=sort(aa,'descend');
        id2 = cl_vn{l1}(:,nc1);
        h1(id2,l1)=ii(1);
    end
end
toc
%%
entropy_prob = zeros(1,N);
if opt_entr==1
    opt1='Entropy';
    for k=1:N/2
        ent0=mne0{1,1}(h1(2*k,1),k);
        ent1=mne1{1,1}(h1(2*k,1),k);
        entropy_prob(2*k-1:2*k)=log2(q)-[ent0 ent1];
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
% save([pwd '\optimized_coefs\N' num2str(N), '_q' num2str(q) '_snr' num2str(SNRs_db+0.01,3) 'dB_observ' num2str(obsrv_cnt) '_' opt1 '.mat'],'h1','sorted_entropy_prob','ch_idx', 'opt1', 'q' )