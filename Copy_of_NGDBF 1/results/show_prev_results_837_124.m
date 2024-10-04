clear

a = load('C:\Users\Abdallah Abdallah\Documents\Abdallah france\OneDrive\matlab-NB-decoders-algorithms\EMS\results\837_124_32_2024_07_24_04_59_08.mat');
b = load('C:\Users\Abdallah Abdallah\Documents\Abdallah france\OneDrive\matlab-NB-decoders-algorithms\AMA MVSF\results\837_124_32_2024_07_23_23_26_04.mat');
c = load('C:\Users\Abdallah Abdallah\Documents\Abdallah france\OneDrive\matlab-NB-decoders-algorithms\NGDBFoctave NB\results\Non_Binary_BPSK_GDBF_837_124_32_dated_2024_07_24_18_12_02.mat');

conf_detail = a.conf_detail;


fields = fieldnames(a.conf_detail);
text_lines = cell(numel(fieldnames(conf_detail)), 1);
for i = 1:numel(fields)
    text_lines{i} = conf_detail.(fields{i});
end

figure;
subplot(1,2,1)
semilogy(a.ebn0, a.FERstat' ,'bo-',b.ebn0, b.FERstat','ro-', c.ebn0, c.FERR','go-', 'Linewidth', 1)
grid on
xlim([0 6])
ylim([1e-6 1e-0])
legend({'FER: EMS nm=16 iter=15';'FER: AMA MVSF 6 attemps of 16 iterations';'FER: NB-GDSF 1000 iters +'}, 'Location', 'southwest');
xlabel('EbNo (dB)')
ylabel('FER')
subplot(1,2,2)

semilogy(a.ebn0, a.BERstat' ,'bo-',b.ebn0, b.BERstat','ro-', c.ebn0, c.BERR','go-','Linewidth', 1)
grid on
xlim([0 6])
ylim([1e-9 1e-0])
legend({'BER: EMS nm=16 iter=15';'BER: AMA MVSF 6 attemps of 16 iterations';'BER: NB-GDSF 1000 iters +'}, 'Location', 'southwest');
xlabel('EbNo (dB)')
ylabel('BER')

