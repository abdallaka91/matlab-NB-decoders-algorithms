clc; clear; close all;

%% Parameters
N = 1e6;        % Number of bits
Rb = 1e3;       % Bit rate (Hz)
B = Rb / 2;     % Signal bandwidth for BPSK
fs = 4*B;       % Sampling frequency (fs >= 2B for Nyquist)
Ts = 1/fs;      % Sampling period
T_symbol = 1/Rb;% Symbol duration
Ns = fs/Rb;     % Samples per symbol

Eb_N0_dB = 0:0.5:5; % SNR values in dB
BER = zeros(1, length(Eb_N0_dB)); % BER storage
CNT = zeros(1, length(Eb_N0_dB)); % BER storage

%% Generate Random Bits
bits = randi([0 1], 1, N);  % Random 0s and 1s

%% BPSK Modulation (Baseband)
tx_signal = repelem(2*bits - 1, Ns);  % BPSK: 0 -> -1, 1 -> +1 (Repeat for oversampling)

%% AWGN Channel
for snr_idx = 1:length(Eb_N0_dB)
    Eb_N0 = 10^(Eb_N0_dB(snr_idx)/10);  % Convert dB to linear scale
    noise_power = (1/(2*Eb_N0)); % Noise power (baseband)
    noise = sqrt(noise_power) * randn(size(tx_signal)); % AWGN

    rx_signal = tx_signal + noise;  % Received signal with noise

    %% Receiver Low-Pass Filtering (Matched Filter)
    h = ones(1, Ns);  % Simple rectangular matched filter
    rx_filtered = conv(rx_signal, h) / Ns;
    rx_filtered = rx_filtered(Ns/2:end-Ns/2);

    %% Symbol Sampling & Decision
    sampled_symbols = rx_filtered(Ns:Ns:end); % Sample at symbol centers
    demod_bits = sampled_symbols > 0; % Decision rule

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
title('BPSK BER Performance over AWGN (Baseband with Filtering)');
legend('Simulated BER');
