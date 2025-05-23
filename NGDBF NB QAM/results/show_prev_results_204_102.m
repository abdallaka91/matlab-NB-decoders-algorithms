clear

a = load('C:\Users\Abdallah Abdallah\Documents\Personals\OneDrive\matlab-NB-decoders-algorithms\EMS_QAM\results\204.102.3.6.16_2024_08_13_02_39_06.mat');
b = load('C:\Users\Abdallah Abdallah\Documents\Personals\OneDrive\matlab-NB-decoders-algorithms\EMS_QAM\results\204.102.3.6.16_2024_08_13_06_05_03.mat');
c = load('C:\Users\Abdallah Abdallah\Documents\Personals\OneDrive\matlab-NB-decoders-algorithms\NGDBF NB QAM\results\Non_Binary_R_QAM_GDBF_204.102.3.6.16_dated_2024_08_13_02_28_30.mat');
d = load('C:\Users\Abdallah Abdallah\Documents\Personals\OneDrive\matlab-NB-decoders-algorithms\NGDBF NB QAM\results\Non_Binary_R_QAM_GDBF_204.102.3.6.16_dated_2024_08_13_02_45_46.mat');
e = load('C:\Users\Abdallah Abdallah\Documents\Personals\OneDrive\matlab-NB-decoders-algorithms\NGDBF NB QAM\results\Non_Binary_R_QAM_GDBF_204.102.3.6.16_dated_2024_08_13_03_54_11.mat');


figure;
subplot(1,2,1)
semilogy(a.Snr_db, a.FERstat' ,'bo-',b.Snr_db, b.FERstat','ro-',...
    c.Snr_db, c.FERR','mo--',d.Snr_db, d.FERR','go--',e.Snr_db, e.FERR','ko--', 'Linewidth', 1)
grid on
xlim([6 14])
ylim([1e-6 1e-0])
legend({'FER: EMS nm=4 iter=15';'FER: EMS nm=8 iter=15';...
    'FER: qam NB-GDBF 1000 iters +';'FER: qam NB-GDBF 2000 iters +' ;...
    'FER: qam NB-GDBF 4000 iters +'}, 'Location', 'southwest');
xlabel('SNR (dB)')
ylabel('FER')

subplot(1,2,2)
semilogy(a.Snr_db, a.BERstat' ,'bo-',b.Snr_db, b.BERstat','ro-',...
    c.Snr_db, c.BERR','mo--',d.Snr_db, d.BERR','go--',e.Snr_db, e.BERR','ko--', 'Linewidth', 1)
grid on
xlim([6 14])
ylim([1e-9 1e-0])
legend({'BER: EMS nm=4 iter=15';'BER: EMS nm=8 iter=15';...
    'BER: qam NB-GDBF 1000 iters +';'BER: qam NB-GDBF 2000 iters +' ;...
    'BER: qam NB-GDBF 4000 iters +'}, 'Location', 'southwest');
xlabel('SNR (dB)')
ylabel('BER')

