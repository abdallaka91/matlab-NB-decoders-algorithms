clear
pth1 = 'C:\Users\Abdallah Abdallah\Documents\Personals\matlab-NB-decoders-algorithms\AMA MVSF\results';

a1=load(fullfile(pth1,'204.102.3.6.16_2024_12_05_15_41_21'));
a2=load(fullfile(pth1,'204.102.3.6.16_2024_12_05_20_21_29'));
a3=load(fullfile(pth1,'204.102.3.6.16_2024_12_05_21_08_43'));
a4=load(fullfile(pth1,'204.102.3.6.16_2024_12_05_22_44_14'));
a5=load(fullfile(pth1,'204.102.3.6.16_2024_06_04_10_49_20'));
figure;

semilogy(a1.ebn0, a1.BERstat)
hold on
semilogy(a2.ebn0, a2.BERstat)
semilogy(a3.ebn0, a3.BERstat)
semilogy(a4.ebn0, a4.BERstat)
hold on
semilogy(a5.ebn0, a5.BERstat)
% semilogy(a6.ebn0, a6.BERstat)
% semilogy(a7.ebn0, a7.BERstat)

hold on

legend({...
    [num2str(a1.max_attempt) ' att, ' num2str(a1.max_iter) ' iter'],...
    [num2str(a2.max_attempt) ' att, ' num2str(a2.max_iter) ' iter'],...
    [num2str(a3.max_attempt) ' att, ' num2str(a3.max_iter) ' iter'],...
    [num2str(a4.max_attempt) ' att, ' num2str(a4.max_iter) ' iter'],...
    [num2str(a5.max_attempt) ' att, ' num2str(a5.max_iter) ' iter']})
grid on
xlim([2 5])
%     
%     [num2str(a5.max_attempt) ' att, ' num2str(a5.max_iter) ' iter'],...
%     [num2str(a6.max_attempt) ' att, ' num2str(a6.max_iter) ' iter'],...

