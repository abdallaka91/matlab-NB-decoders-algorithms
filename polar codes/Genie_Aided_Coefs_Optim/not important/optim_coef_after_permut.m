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

SNRs_db = -5;
SNRs = 10.^(SNRs_db/10);
N0 = 1./(SNRs);
sigma = sqrt(N0);
words = (0:q-1);
alph_bin =  fliplr(dec2bin(words, p) - 48);
alph_bin_mod = (-1).^alph_bin;
h01=1;

sigm = sigma(1);
opt_entr=1;
n_obs=500;

perms_n=cell(n,1);

for i=n:-1:1
 perms_n{i}=nan(2^(n-i),q);
end
perms_n_inv = perms_n;

mne0=cell(n,1);
mne1=mne0;
e0=mne0;
e1=mne1;
Deltas=mne0;
tic

% for i=1:n-1
%  coefs{i}=64+2^i+(0:2^(i-1)-1);
% end

% 0 dB
% coefsp={[];...
%  [54 46];...
%  [24 60 51 30];...
%  [15 61 29 21  8 56 56 57];...
%  [43 23 61 44 42 13  8 16 41 19 60  7 37  7 37 47];...
%  [16 19 36 15 6 6 21 60 60 17 13 28 4 3 24 46 3 16 60 6 6 58 25 63 45 13 27 20 24 61 63 1];...
%  [15 27 3 46 48 41 50 32 43 7 16 12 10 57 14 32 25 16 41 59 16 42 32 42 42 32 32 58 50 13 53 48 50 32 36 9 52 43 43 16 28 60 31 58 9 46 9 58 56 22 23 28 43 16 37 1 44 46 53 1 49 1 1 1]
%  [29 18 42 51 2 26 60 53 38 52 11 15 55 53 57 4 13 30 36 25 59 38 35 59 51 37 16 26 14 24 15 57 3 55 27 38 54 1 4 1 62 34 5 5 19 1 24 25 57 16 16 27 23 5 1 47 5 48 60 47 47 7 18 46 48 55 31 5 53 5 45 54 5 49 30 29 14 43 30 29 6 52 2 10 60 56 19 56 29 43 53 19 43 8 30 1 11 26 37 55 15 49 3 1 50 45 26 1 22 1 1 1 10 19 34 1 22 1 1 1 10 1 1 1 1 1 1 1]}; % previous layer optimized coefs

%5 dB
% coefsp={[];...
%     [                                                                            32 1];...
%     [                                                                     36 48 48 23];...
%     [                                                            9 33 25 25 25 36 5 1];...
%     [                                       20 48 13 34 21 59 41 11 7 54 61 1 3 1 1 1];...
%     [60 32 58 52 53 12 48 41 59 44 57 53 27 29 16 1 2 10 63 31 46 1 1 1 25 1 1 1 1 1 1 1];...
%     };
coefsp={[]};
nnn=length(coefsp);

% uu=load('C:\Users\Abdallah Abdallah\Documents\Personals\matlab-NB-decoders-algorithms\polar codes\Genie_Aided_Coefs_Optim\optimized_coefs\N64_q64_snr5dB_observ500_Entropy_all_param.mat');
% coefsp=uu.coefs;

perms_n{n}=[56 1 20 55 31 43 2 44 51 29 46 26 37 6 9 48 45 14 27 34 54 24 49 5 0 52 7 41 42 19 28 63 10 36 13 57 17 50 39 30 22 47 32 3 60 8 59 21 35 23 62 16 4 61 40 11 25 58 53 12 15 33 18 38];
perms_n_inv{n}(perms_n{n}+1)=0:q-1;
coefs=cell(n,1);

coefs(1:nnn)=coefsp(1:nnn);

for n0 = 2:nnn
 N1=2^n0;
 for k=1:N1/2
  perms_n{n-n0+1}(k,:)=mul_mat(1+coefs{n0}(k),:);
  perms_n_inv{n-n0+1}(k,perms_n{n-n0+1}(k,:)+1)=0:q-1;
 end

end
for n0 = nnn+1:n
 perms_n0=perms_n(n-n0+1:n);
 perms_n0_inv=perms_n_inv(n-n0+1:n);
 N1=2^n0;
 mne0{n-n0+1}=nan(q-1,2^(n0-1));
 mne1{n-n0+1}=nan(q-1,2^(n0-1));
 e0{n-n0+1}=nan(q-1,n_obs, 2^(n0-1));
 e1{n-n0+1}=nan(q-1,n_obs, 2^(n0-1));


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

 h1 = nan(N1,n0);
 for i=2:n0-1
  for j=1:2^(i-1)
   h1(cl_vn{n0-i+1}(:,j),n0-i+1)=coefs{i}(j);
  end
 end

 for i = 1 : n_obs
  PP=nan(q,N1, n0+1);
  info_seq=randi([0 q-1], 1, N1);
  %   info_seq=[10 38 3 38];
  info_seq=info_seq(1:N1);
  q1=q-1;
  h0=0;
  while h0<q-1
   h0=h0+1;
   if n0>1
    for k=1:N1/2
     perms_n0{1}(k,:)=mul_mat(1+h0,:);
     perms_n0_inv{1}(k,perms_n0{1}(k,:)+1)=0:q-1;
     %      perms_n0{1}(k,:)=uuu;
    end

   end
   [info_seq,x,m1, y]=gen_noisy_seq1(N1,q,sigm,perms_n0, add_mat, alph_bin_mod, 1, info_seq);
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
       [Pr0, Pr1]= F1_permutation(L1, L2, add_mat,perms_n0{l2}(nc1,:), u0);
       PP(:,ii1,l2)=Pr0;
       PP(:,ii2,l2)=Pr1;
      else
       L1=squeeze(PP(:,ii1,1+1));
       L2=squeeze(PP(:,ii2,1+1));
       [Pr0, Pr1]= F1_permutation(L1, L2, add_mat,perms_n0{l2}(nc1,:), u0);
       %        if n0==1 && nc1==1
       %        ttt=ttt+Pr1;
       %        tttt=tttt+L2;
       %        end
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
  aa=Deltas{n-n0+1}(:, nc1);
  [~,ii]=sort(aa,'descend');
  if n0>1
   coefs{n0}(nc1)=ii(1);
   perms_n{n-n0+1}(nc1,:)=mul_mat(ii(1)+1,:);
   h1(cl_vn{1}(nc1), 1)=ii(1);
  end
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
  entropy_prob(2*k-1:2*k) = [ent0 ent1];
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

for l1=2:-1:1
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

% save([pwd '\optimized_coefs\N' num2str(N), '_q' num2str(q) '_snr' num2str(SNRs_db+0.01,3) 'dB_observ' num2str(n_obs) '_' opt1 '.mat'],'h1','sorted_entropy_prob','ch_idx', 'opt1', 'q' )