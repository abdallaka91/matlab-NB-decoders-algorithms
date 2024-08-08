% clear
% pth6 = (fullfile(pwd, 'results/'));
% a = load(fullfile(pth6, '204.102.3.6.16_2024_05_30_14_21_41.mat'));
% b = load(fullfile(pth6, '204.102.3.6.16_2024_05_29_17_49_10.mat'));
% c = load(fullfile(pth6, '204.102.3.6.16_2024_05_29_17_49_10.mat'));
% 
% %%
% conf_detail = a.conf_detail;
% 
% 
% fields = fieldnames(a.conf_detail);
% text_lines = cell(numel(fieldnames(conf_detail)), 1);
% for i = 1:numel(fields)
%     text_lines{i} = conf_detail.(fields{i});
% end
% 
% ebn00 = 1.5:0.25:3;
% fer00 = [6e-1 2.5e-1 1e-1 2.8e-2 3.2e-3 3e-4 1.8e-5];
% ber00 = [7e-2 3e-2 1e-2 2.6e-3 3e-4 3e-5 2e-6];
% 
% figure;
% semilogy(a.ebn0, a.FERstat' ,'o-',b.ebn0, b.FERstat','o-')
% hold on
% semilogy(ebn00, fer00 ,'*-')
% grid on
% xlim([0 6])
% ylim([1e-7 1e-0])
% % legend({'12 iters - 1 attempt', '96 iters - 1 attempt', '12 iters - 8 attempts'})
% legend({'MV-SF | 20 iters | 1 attempt', 'MV-SF | 20 iters | 4 attempts', 'EMS (from the paper)'})
% xlabel('EbNo (dB)')
% ylabel('FER')
% % position = [0.1 + 0.4 * 8 / 100, 0.8 - 0.8 * 50 / 100, 0.8, 0.1];
% % annotation('textbox', position, 'String', text_lines, 'Interpreter', 'none', 'FitBoxToText', 'on', 'EdgeColor', 'none');
% % fig_pos = get(gcf, 'Position');
% % fig_pos(3:4) = fig_pos(3:4) .* [1.4, 1.2];
% % set(gcf, 'Position', fig_pos);
% %%
% ebn00 = 1.5:0.25:3;
% fer00 = [6e-1 2.5e-1 1e-1 2.8e-2 3.2e-3 3e-4 1.8e-5];
% ber00 = [7e-2 3e-2 1e-2 2.6e-3 3e-4 3e-5 2e-6];



clear
pth6 = (fullfile(pwd, 'results/'));
a = load(fullfile(pth6, '204.102.3.6.16_2024_05_30_14_21_41.mat'));
b = load(fullfile(pth6, '204.102.3.6.16_2024_05_29_17_49_10.mat'));
c = load(fullfile(pth6, '204.102.3.6.16_2024_05_29_17_49_10.mat'));

%%
conf_detail = a.conf_detail;


fields = fieldnames(a.conf_detail);
text_lines = cell(numel(fieldnames(conf_detail)), 1);
for i = 1:numel(fields)
    text_lines{i} = conf_detail.(fields{i});
end

ebn00 = 1.5:0.25:3;
fer00 = [6e-1 2.5e-1 1e-1 2.8e-2 3.2e-3 3e-4 1.8e-5];
ber00 = [7e-2 3e-2 1e-2 2.6e-3 3e-4 3e-5 2e-6];

figure;
semilogy(a.ebn0, a.BERstat' ,'o-',b.ebn0, b.BERstat','o-')
hold on
semilogy(ebn00, ber00 ,'*-')
grid on
xlim([0 6])
ylim([1e-7 1e-0])
% legend({'12 iters - 1 attempt', '96 iters - 1 attempt', '12 iters - 8 attempts'})
legend({'MV-SF | 20 iters | 1 attempt', 'MV-SF | 20 iters | 4 attempts', 'EMS (from the paper)'})
xlabel('EbNo (dB)')
ylabel('BER')
% position = [0.1 + 0.4 * 8 / 100, 0.8 - 0.8 * 50 / 100, 0.8, 0.1];
% annotation('textbox', position, 'String', text_lines, 'Interpreter', 'none', 'FitBoxToText', 'on', 'EdgeColor', 'none');
% fig_pos = get(gcf, 'Position');
% fig_pos(3:4) = fig_pos(3:4) .* [1.4, 1.2];
% set(gcf, 'Position', fig_pos);
%%
ebn00 = 1.5:0.25:3;
fer00 = [6e-1 2.5e-1 1e-1 2.8e-2 3.2e-3 3e-4 1.8e-5];
ber00 = [7e-2 3e-2 1e-2 2.6e-3 3e-4 3e-5 2e-6];