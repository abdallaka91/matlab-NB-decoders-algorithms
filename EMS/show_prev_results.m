clear
pth6 = (fullfile(pwd, 'results/'));
a = load(fullfile(pth6, '204.102.3.6.16_2024_06_04_10_49_20_MVSF.mat'));
b = load(fullfile(pth6, '204.102.3.6.16_2024_07_04_04_53_06.mat'));
c = load(fullfile(pth6, '204.102.3.6.16_2024_07_04_23_35_28.mat'));

conf_detail = a.conf_detail;


fields = fieldnames(a.conf_detail);
text_lines = cell(numel(fieldnames(conf_detail)), 1);
for i = 1:numel(fields)
    text_lines{i} = conf_detail.(fields{i});
end

figure;
semilogy(a.ebn0, a.FERstat' ,'o-',b.ebn0, b.FERstat','o-', c.ebn0, c.FERstat','o-')
grid on
xlim([0 6])
ylim([1e-7 1e-0])
legend({'BeiDou 88-44-GF64', '204-102-GF16', '255-176-GF256'})
xlabel('EbNo (dB)')
ylabel('FER')

