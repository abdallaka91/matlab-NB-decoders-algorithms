clear

a = load('C:\Users\Abdallah Abdallah\Documents\Personals\OneDrive\matlab-NB-decoders-algorithms\EMS\results\204.102.3.6.16_2024_07_24_21_24_18 - Copy.mat');
b = load('C:\Users\Abdallah Abdallah\Documents\Personals\OneDrive\matlab-NB-decoders-algorithms\AMA MVSF\results\204.102.3.6.16_2024_06_04_10_49_20.mat');
c = load('C:\Users\Abdallah Abdallah\Documents\Personals\OneDrive\matlab-NB-decoders-algorithms\AMA MVSF\results\204.102.3.6.16_2024_06_03_16_52_09.mat');
d = load('C:\Users\Abdallah Abdallah\Documents\Personals\OneDrive\matlab-NB-decoders-algorithms\NGDBF NB\results\Non_Binary_BPSK_GDBF_204.102.3.6.16_dated_2024_07_24_22_55_10 - Copy.mat');


figure;
subplot(1,2,1)
semilogy(a.ebn0, a.FERstat' ,'bo-',b.ebn0, b.FERstat','ro-', c.ebn0, c.FERstat','mo-', d.ebn0, d.FERR','go-', 'Linewidth', 1)
grid on
xlim([0 6])
ylim([1e-6 1e-0])
legend({'FER: EMS nm=16 iter=15';'FER: AMA MVSF 6 attemps of 16 iterations';...
    'FER: MVSF 96 iterations'; 'FER: NB-GDSF 1000 iters +'}, 'Location', 'southwest');
xlabel('EbNo (dB)')
ylabel('FER')

subplot(1,2,2)
semilogy(a.ebn0, a.BERstat' ,'bo-',b.ebn0, b.BERstat','ro-', c.ebn0, c.BERstat','mo-', d.ebn0, d.BERR','go-','Linewidth', 1)
grid on
xlim([0 6])
ylim([1e-9 1e-0])
legend({'BER: EMS nm=16 iter=15';'BER: AMA MVSF 6 attemps of 16 iterations';...
    'BER: MVSF 96 iterations';'BER: NB-GDSF 1000 iters +'}, 'Location', 'southwest');
xlabel('EbNo (dB)')
ylabel('BER')

