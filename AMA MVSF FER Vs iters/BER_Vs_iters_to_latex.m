clear
pth = 'C:\Users\Abdallah Abdallah\Documents\Personals\matlab-NB-decoders-algorithms\AMA MVSF FER Vs iters\results';


a = load(fullfile(pth, '204.102.3.6.16_2024_12_11_14_44_00.mat'));
b = load(fullfile(pth, '204.102.3.6.16_2024_12_11_14_27_21.mat'));
% a = load(fullfile(pth, 'BeiDou_44_88_GF64_2024_12_11_14_26_39.mat'));
% b = load(fullfile(pth, 'BeiDou_44_88_GF64_2024_12_11_14_45_41.mat'));
fl_nm1 = ['C:\Users\Abdallah Abdallah\Documents\Personals\OneDrive\shared PhD\paper drafts\6672d05ce10082c5b1416aa3\BER\' a.H_matrix_mat_fl_nm '_BER_vs_iters_attempts.txt'];

% b.BER_iter1(15)=b.BER_iter1(15)*1.05;
% b.BER_iter1(16)=b.BER_iter1(16)*1.05;

outputTable = table((1:length(b.BER_iter1))', (a.BER_iter1)', (a.FER_iter1)', (b.BER_iter1)',(b.FER_iter1)',...
    'VariableNames', {'iters', 'multiple_att_ber','multiple_att_fer', 'single_att_ber', 'single_att_fer'});

% Define output file name

% Write table to file with specified format
writetable(outputTable, fl_nm1, 'Delimiter', '\t', 'WriteVariableNames', true);

figure;
semilogy(a.BER_iter1)
hold on

semilogy(b.BER_iter1)
grid on