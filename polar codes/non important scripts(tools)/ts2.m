clear

% rng(0)
filePath = 'C:\Users\Abdallah Abdallah\Documents\Personals\matlab-NB-decoders-algorithms\polar codes\december 2024\bpsk_ccsk_gf64\biawgn_ccsk[6]_gfdim6_n4.plr';
[p, data] = read_nb_polar_files_parameter(filePath);
% best channels to worst channels indexes
% load('C:\Users\Abdallah Abdallah\Documents\Personals\matlab-NB-decoders-algorithms\polar codes\december 2024\coefs_G_Ginv_CCSK\biawgn_ccsk[6]_gfdim6_n4.mat');

i2=25;
p=6;
q=2^p;
mainPath = pwd;
% mainPath = fileparts(projectPath);
related_variables_pth = fullfile(mainPath, 'related_variables');
pth3 = (fullfile(related_variables_pth, 'GF_arithm'));
fl_nm = ['arith_' num2str(q) '.mat'];
load(fullfile(pth3, fl_nm));

N=length(data(1).Pe);

coefs = data(i2).gf_coef;
u=randi([0 q-1],1,N)';
u=(0:N-1)';
[coefs1,I1] = coef_2_coef(coefs);
 [x1,x2,I1]=encode2(u, coefs,add_mat, mul_mat);
% % y3(I1+1,:)=y2;
% x = encode1(u, coefs,add_mat, mul_mat);
%  u1 = decode1(x, coefs,add_mat, div_mat);