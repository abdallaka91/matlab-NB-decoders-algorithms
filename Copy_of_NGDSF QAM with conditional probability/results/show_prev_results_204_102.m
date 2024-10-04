clear

a = load('/home/matlab-runner/abdallah_files/EMS_QAM/results/204.102.3.6.16_2024_08_13_02_39_06.mat');
b = load('/home/matlab-runner/abdallah_files/EMS_QAM/results/204.102.3.6.16_2024_08_13_06_05_03.mat');
c = load('/home/matlab-runner/abdallah_files/NGDBF NB QAM/results/Non_Binary_R_QAM_GDBF_204.102.3.6.16_dated_2024_08_13_02_28_30.mat');
d = load('/home/matlab-runner/abdallah_files/NGDBF NB QAM/results/Non_Binary_R_QAM_GDBF_204.102.3.6.16_dated_2024_08_13_02_45_46.mat');
e = load('/home/matlab-runner/abdallah_files/NGDBF NB QAM/results/Non_Binary_R_QAM_GDBF_204.102.3.6.16_dated_2024_08_13_03_54_11.mat');
f = load('/home/matlab-runner/abdallah_files/NGDSF QAM with conditional probability/results/Non_Binary_RQAM_GD_SymboleFlipping_204.102.3.6.16_dated_2024_09_16_00_51_49.mat');
g = load('/home/matlab-runner/abdallah_files/NGDSF QAM with conditional probability/results/Non_Binary_RQAM_GD_SymboleFlipping_204.102.3.6.16_dated_2024_09_16_01_09_36.mat');
h = load('/home/matlab-runner/abdallah_files/NGDSF QAM with conditional probability/results/Non_Binary_RQAM_GD_SymboleFlipping_204.102.3.6.16_dated_2024_09_16_23_19_14.mat');



lnw = 2;

figure;
subplot(1,2,1)
semilogy(a.Snr_db, a.BERstat' ,'b*-', 'Linewidth', lnw);
hold on
semilogy(b.Snr_db, b.FERstat','r*-', 'Linewidth', lnw);
semilogy(c.Snr_db, c.FERR','mo--', 'Linewidth', lnw);
semilogy(d.Snr_db, d.FERR','go--', 'Linewidth', lnw);
semilogy(e.Snr_db, e.FERR','co--', 'Linewidth', lnw);
semilogy(f.Snr_db, f.FERR','k+--', 'Linewidth', lnw);
semilogy(g.Snr_db, g.FERR','+--', 'Linewidth', lnw);
semilogy(h.Snr_db, h.FERR','+--', 'Linewidth', lnw);


grid on
xlim([6 13])
ylim([1e-7 1e-0])
legend({'FER: EMS nm=4 iter=15';'FER: EMS nm=8 iter=15';...
    'FER: qam NB-GDBF 1000 iters +';'FER: qam NB-GDBF 2000 iters +' ;...
    'FER: qam NB-GDBF 300 iters +';...
    'FER: qam NB-GD-SymbFlip 1000 iters +';...
    'FER: qam NB-GD-SymbFlip 2000 iters +';...
    'FER: qam NB-GD-SymbFlip 4000 iters +'}, 'Location', 'southwest');
xlabel('SNR (dB)')
ylabel('FER')
subplot(1,2,2)
semilogy(a.Snr_db, a.BERstat' ,'b*-', 'Linewidth', lnw);
hold on
semilogy(b.Snr_db, b.BERstat','r*-', 'Linewidth', lnw);
semilogy(c.Snr_db, c.BERR','mo--', 'Linewidth', lnw);
semilogy(d.Snr_db, d.BERR','go--', 'Linewidth', lnw);
semilogy(e.Snr_db, e.BERR','co--', 'Linewidth', lnw);
semilogy(f.Snr_db, f.BERR','k+--', 'Linewidth', lnw);
semilogy(g.Snr_db, g.BERR','+--', 'Linewidth', lnw);
semilogy(h.Snr_db, h.BERR','+--', 'Linewidth', lnw);


grid on
xlim([6 13])
ylim([1e-8 1e-0])
legend({'BER: EMS nm=4 iter=15';'BER: EMS nm=8 iter=15';...
    'BER: qam NB-GDBF 1000 iters +';'BER: qam NB-GDBF 2000 iters +' ;...
    'BER: qam NB-GDBF 4000 iters +';...
    'BER: qam NB-GD-SymbFlip 1000 iters +';...
    'BER: qam NB-GD-SymbFlip 2000 iters +';...
    'BER: qam NB-GD-SymbFlip 300 iters +'}, 'Location', 'southwest');
xlabel('SNR (dB)')
ylabel('BER')