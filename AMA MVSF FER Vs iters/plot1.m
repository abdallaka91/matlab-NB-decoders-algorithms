clear
pth = 'C:\Users\Abdallah Abdallah\Documents\Personals\matlab-NB-decoders-algorithms\AMA MVSF FER Vs iters\results';

a = load(fullfile(pth, '204.102.3.6.16_2024_12_02_11_28_52.mat'));
b = load(fullfile(pth, '204.102.3.6.16_2024_12_02_11_34_52.mat'));
c = load(fullfile(pth, '204.102.3.6.16_2024_12_02_12_40_38.mat'));

figure(3)
semilogy(a.FER_iter1)
hold on
semilogy(b.FER_iter1, 'o-')
% semilogy(c.FER_iter1, '*-')

legend({[num2str(a.max_attempt) ' att, ' num2str(a.max_iter) ' iters'],...
    [num2str(b.max_attempt) ' att, ' num2str(b.max_iter) ' iters']})%,...
%     [num2str(c.max_attempt) ' att, ' num2str(c.max_iter) ' iters']})


figure(4)
semilogy(a.BER_iter1)
hold on
semilogy(b.BER_iter1, 'o-')
% semilogy(c.BER_iter1, '*-')

legend({[num2str(a.max_attempt) ' att, ' num2str(a.max_iter) ' iters'],...
    [num2str(b.max_attempt) ' att, ' num2str(b.max_iter) ' iters']})%,...
%     [num2str(c.max_attempt) ' att, ' num2str(c.max_iter) ' iters']})