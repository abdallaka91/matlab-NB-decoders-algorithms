clear

% rng(0)
filePath = 'C:\Users\Abdallah Abdallah\Documents\Personals\matlab-NB-decoders-algorithms\polar codes\december 2024 - Copy\bpsk_gf64\biawgn_gfdim6_n8.plr';
[p, data] = read_nb_polar_files_parameter(filePath);
% best channels to worst channels indexes
q=2^p;
mainPath = pwd;
% mainPath = fileparts(projectPath);
related_variables_pth = fullfile(mainPath, 'related_variables');
pth3 = (fullfile(related_variables_pth, 'GF_arithm'));
fl_nm = ['arith_' num2str(q) '.mat'];
load(fullfile(pth3, fl_nm));

N=length(data(1).Pe);
    
[~, fileName, ~] = fileparts(filePath);

folder_of_mat_file = fullfile(pwd, 'data_bpsk');

save(fullfile(folder_of_mat_file, [fileName '.mat']), 'data')



