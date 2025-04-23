clc; clear; close all;

%% Parameters
N = 1e4;        % Number of bits
Rb = 1e3;       % Bit rate (Hz)
B = Rb / 2;     % Signal bandwidth for BPSK
fc = 5e3;       % Carrier frequency (Hz)
fs = 10*fc;     % Sampling frequency (fs >> 2B, practical case)
Ts = 1/fs;      % Sampling period
T_symbol = 1/Rb;% Symbol duration
Ns = fs/Rb;     % Samples per symbol

Eb_N0_dB = -20:0.25:15; % SNR values in dB
BER = zeros(1, length(Eb_N0_dB)); % BER storage
CNT = zeros(1, length(Eb_N0_dB)); % Error count storage

%% Generate Random Bits
bits = randi([0 1], 1, N);  % Random 0s and 1s

%% BPSK Baseband Modulation
baseband_signal = repelem(2*bits - 1, Ns);  % BPSK: 0 -> -1, 1 -> +1

%% Carrier Modulation (RF)
t = (0:length(baseband_signal)-1) * Ts; % Time vector
carrier = cos(2*pi*fc*t);   % RF carrier
tx_signal = baseband_signal .* carrier;  % Modulated RF signal

%% AWGN Channel (Noise added to RF Signal)
for snr_idx = 1:length(Eb_N0_dB)
    Eb_N0 = 10^(Eb_N0_dB(snr_idx)/10);  % Convert dB to linear scale
    noise_power = (1/(2*Eb_N0)); % Noise power
    noise = sqrt(noise_power) * randn(size(tx_signal)); % AWGN

    rx_signal = tx_signal + noise;  % Received signal with noise

    %% Coherent Demodulation (Multiplication with Local Carrier)
    rx_mixed = rx_signal .* (2 * cos(2*pi*fc*t)); % Multiply by 2cos(wt)

    %% Low-Pass Filtering (Matched Filter)
    h = ones(1, Ns/2);  % Simple rectangular filter
    rx_filtered = cconv(rx_mixed, h, length(h)+length(rx_mixed)-1) / Ns;
    rx_filtered = rx_filtered(length(h)/2:end-length(h)/2);

    %% Symbol Sampling & Decision
    sampled_symbols = rx_filtered(Ns/2:Ns:end); % Sample at symbol centers
    demod_bits = sampled_symbols > 0; % Decision rule: Positive -> 1, Negative -> 0

    %% BER Calculation
    CNT(snr_idx) = sum(bits ~= demod_bits);
    BER(snr_idx) = CNT(snr_idx)  / N;
end

%% Plot BER vs Eb/N0
figure;
semilogy(Eb_N0_dB, BER, 'bo-', 'LineWidth', 2);
grid on;
xlabel('E_b/N_0 (dB)');
ylabel('Bit Error Rate (BER)');
title('BPSK BER Performance over AWGN (RF Modulation & Coherent Demodulation)');
legend('Simulated BER');
