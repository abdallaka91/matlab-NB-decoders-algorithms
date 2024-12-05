clear
fl_nm = 'C:\Users\Abdallah Abdallah\Documents\Personals\matlab-NB-decoders-algorithms\AMA MVSF\results\204.102.3.6.16_2024_10_29_15_54_13.mat';
fl_nm1 = 'C:\Users\Abdallah Abdallah\Documents\Personals\OneDrive\shared PhD\paper drafts\AMA MVSF bib\204_102_12_deviations_v3_etta_3_16iters_1_attmps.txt';

data = load(fl_nm);  % Replace 'yourfile.mat' with the actual filename

NN1 = 1;
NN2 = 12;
Ebn0 = data.ebn0(NN1:NN2)';      % Signal-to-noise ratio (SNR)
BERstat = data.BERstat(NN1:NN2); % Bit Error Rate (BER)
FERstat = data.FERstat(NN1:NN2); % Frame Error Rate (FER)

% Create table with headers
outputTable = table(Ebn0, BERstat, FERstat, 'VariableNames', {'SNR', 'BER', 'FER'});

% Define output file name

% Write table to file with specified format
% writetable(outputTable, fl_nm1, 'Delimiter', '\t', 'WriteVariableNames', true);

% Display message
disp(['Data written to ', fl_nm1]);


figure(1)
semilogy(Ebn0, BERstat)
hold on
grid on