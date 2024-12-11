clear
pth = 'C:\Users\Abdallah Abdallah\Documents\Personals\matlab-NB-decoders-algorithms\AMA MVSF FER Vs iters\results';

a = load(fullfile(pth, 'BeiDou_44_88_GF64_2024_12_03_14_43_23.mat'));
b = load(fullfile(pth, 'BeiDou_44_88_GF64_2024_12_03_14_50_29.mat'));
c = load(fullfile(pth, 'BeiDou_44_88_GF64_2024_12_03_14_51_44.mat'));

figure(1)
semilogy(a.FER_iter1)
hold on
semilogy(b.FER_iter1, 'o-')
% semilogy(c.FER_iter1, '*-')

legend({[num2str(a.max_attempt) ' att, ' num2str(a.max_iter) ' iters'],...
    [num2str(b.max_attempt) ' att, ' num2str(b.max_iter) ' iters']})%,...
%     [num2str(c.max_attempt) ' att, ' num2str(c.max_iter) ' iters']})
xlabel('')


figure(2)
semilogy(a.BER_iter1)
hold on
semilogy(b.BER_iter1, 'o-')
% semilogy(c.BER_iter1, '*-')

legend({[num2str(a.max_attempt) ' att, ' num2str(a.max_iter) ' iters'],...
    [num2str(b.max_attempt) ' att, ' num2str(b.max_iter) ' iters']})%,...
%     [num2str(c.max_attempt) ' att, ' num2str(c.max_iter) ' iters']})