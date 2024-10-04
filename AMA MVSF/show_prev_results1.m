clear
pth6 = (fullfile(pwd, 'results/'));
a = load(fullfile(pth6, 'BeiDou_44_bb_GF64_2024_07_12_05_37_53.mat'));
b = load(fullfile(pth6, 'BeiDou_44_bb_GF64_2024_07_12_03_39_10.mat'));
c = load(fullfile(pth6, 'BeiDou_44_bb_GF64_2024_07_12_00_01_28.mat'));
% d = load(fullfile(pth6, '204.102.3.6.16_2024_06_04_10_49_20.mat'));

%%
conf_detail = a.conf_detail;


fields = fieldnames(a.conf_detail);
text_lines = cell(numel(fieldnames(conf_detail)), 1);
for i = 1:numel(fields)
    text_lines{i} = conf_detail.(fields{i});
end


figure;
semilogy(a.ebn0, a.FERstat' ,'o-')
hold on
semilogy(b.ebn0, b.FERstat','o-')
semilogy(c.ebn0, c.FERstat','o-')
% semilogy(d.ebn0, d.FERstat','o-')
grid on
xlim([0 6])
ylim([1e-7 1e-0])
legend({'FER without perturbation', 'FER with perturbation'})
position = [0.1 + 0.4 * 8 / 100, 0.8 - 0.8 * 50 / 100, 0.8, 0.1];
annotation('textbox', position, 'String', text_lines, 'Interpreter', 'none', 'FitBoxToText', 'on', 'EdgeColor', 'none');
fig_pos = get(gcf, 'Position');
fig_pos(3:4) = fig_pos(3:4) .* [1.4, 1.2];
set(gcf, 'Position', fig_pos);