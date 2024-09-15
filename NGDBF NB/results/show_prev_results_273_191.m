clear

a = load('C:\Users\Abdallah Abdallah\Documents\Personals\OneDrive\matlab-NB-decoders-algorithms\EMS\results\273.191.4.14.64_2024_07_25_23_51_07.mat');
b = load('C:\Users\Abdallah Abdallah\Documents\Personals\OneDrive\matlab-NB-decoders-algorithms\EMS\results\273.191.4.14.64_2024_07_26_04_27_48.mat');
c = load('C:\Users\Abdallah Abdallah\Documents\Personals\OneDrive\matlab-NB-decoders-algorithms\NGDBF NB\results\Non_Binary_BPSK_GDBF_273.191.4.14.64_dated_2024_07_25_15_40_47.mat');
d = load('C:\Users\Abdallah Abdallah\Documents\Personals\OneDrive\matlab-NB-decoders-algorithms\NGDBF NB\results\Non_Binary_BPSK_GDBF_273.191.4.14.64_dated_2024_07_25_23_57_58.mat');
e = load('C:\Users\Abdallah Abdallah\Documents\Personals\OneDrive\matlab-NB-decoders-algorithms\NGDBF NB\results\Non_Binary_BPSK_GDBF_273.191.4.14.64_dated_2024_07_26_04_33_28.mat');
f = load('C:\Users\Abdallah Abdallah\Documents\Personals\OneDrive\matlab-NB-decoders-algorithms\AMA MVSF\results\273.191.4.14.64_2024_07_26_05_37_00.mat');


figure;
subplot(1,2,1)
semilogy(a.ebn0, a.FERstat' ,'ro-',b.ebn0, b.FERstat','bo-', 'Linewidth', 1)
hold on
semilogy(c.ebn0, c.FERR' ,'r*-',d.ebn0, d.FERR','b*-',e.ebn0, e.FERR','g*-', 'Linewidth', 1)
semilogy(f.ebn0, f.FERstat' ,'k+--', 'Linewidth', 1)

grid on
xlim([1 5])
ylim([1e-6 1e-0])
legend({'FER: EMS nm=16 iter=50';'FER: EMS nm=16 iter=15';...
    'FER: NB-GDSF 1000 iters +';'FER: NB-GDSF 2000 iters +';'FER: NB-GDSF 4000 iters +'; ...
    'FER: MVSF iter=50'},...
    'Location', 'southwest');
xlabel('EbNo (dB)')
ylabel('FER')


subplot(1,2,2)

semilogy(a.ebn0, a.BERstat' ,'ro-',b.ebn0, b.BERstat','bo-', 'Linewidth', 1)
hold on
semilogy(c.ebn0, c.BERR' ,'r*-',d.ebn0, d.BERR','b*-',e.ebn0, e.BERR','g*-', 'Linewidth', 1)
semilogy(f.ebn0, f.BERstat' ,'k+--', 'Linewidth', 1)
grid on
xlim([1 5])
ylim([1e-8 1e-1])
legend({'BER: EMS nm=16 iter=50';'BER: EMS nm=16 iter=15';...
    'BER: NB-GDSF 1000 iters +';'BER: NB-GDSF 2000 iters +';'BER: NB-GDSF 4000 iters +'; ...
    'BER: MVSF iter=50'},...
    'Location', 'southwest');
xlabel('EbNo (dB)')
ylabel('BER')

