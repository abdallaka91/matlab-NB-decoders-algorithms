clear
pth1 = 'C:\Users\Abdallah Abdallah\Documents\Personals\matlab-NB-decoders-algorithms\AMA MVSF\results';

a1=load(fullfile(pth1,'BeiDou_44_88_GF64_2024_12_05_15_34_21'));
a2=load(fullfile(pth1,'BeiDou_44_88_GF64_2024_12_05_16_08_55'));
a3=load(fullfile(pth1,'BeiDou_44_88_GF64_2024_12_05_20_21_08'));
a4=load(fullfile(pth1,'BeiDou_44_88_GF64_2024_12_05_20_50_10'));
a5=load(fullfile(pth1,'BeiDou_44_88_GF64_2024_12_05_21_02_25'));
a6=load(fullfile(pth1,'BeiDou_44_88_GF64_2024_12_05_21_16_42'));
a7=load(fullfile(pth1,'BeiDou_44_88_GF64_2024_11_04_18_38_47'));

figure;

semilogy(a1.ebn0, a1.BERstat)
hold on
semilogy(a2.ebn0, a2.BERstat)
semilogy(a3.ebn0, a3.BERstat)
semilogy(a4.ebn0, a4.BERstat)
hold on
semilogy(a5.ebn0, a5.BERstat)
semilogy(a6.ebn0, a6.BERstat)
semilogy(a7.ebn0, a7.BERstat)

hold on

legend({[num2str(a1.max_attempt) ' att, ' num2str(a1.max_iter) ' iter'],...
    [num2str(a2.max_attempt) ' att, ' num2str(a2.max_iter) ' iter'],...
    [num2str(a3.max_attempt) ' att, ' num2str(a3.max_iter) ' iter'],...
    [num2str(a4.max_attempt) ' att, ' num2str(a4.max_iter) ' iter'],...
    [num2str(a5.max_attempt) ' att, ' num2str(a5.max_iter) ' iter'],...
    [num2str(a6.max_attempt) ' att, ' num2str(a6.max_iter) ' iter'],...
    [num2str(a7.max_attempt) ' att, ' num2str(a7.max_iter) ' iter']})
grid on
xlim([2 5])
%%
plot(a1.ebn0,a1.aver_iter)
hold on
plot(a4.ebn0,a4.aver_iter)
grid on


