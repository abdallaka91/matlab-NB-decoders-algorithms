clc; clear; close all;

%% Parameters
N = 1e6;        % Number of bits
Rb = 1e3;       % Bit rate (Hz)
% B = Rb / 2;     % Signal bandwidth for QPSK
% fs = 4*B;       % Sampling frequency (fs >= 2B for Nyquist)

Ts = 1/fs;      % Sampling period
T_symbol = 1/Rb;% Symbol duration
Ns = fs/Rb;     % Samples per symbol

Eb_N0_dB = 0:0.5:5; % SNR values in dB
BER = zeros(1, length(Eb_N0_dB)); % BER storage
CNT = zeros(1, length(Eb_N0_dB)); % BER storage

%% Generate Random Bits for 16-QPSK (4 bits per symbol)
bits = randi([0 1], 1, N);  % Random 0s and 1s

% Group bits into 4-bit symbols for 16-QPSK
symbols = reshape(bits, 4, [])';  % Each row is a 4-bit symbol

% Map 4 bits to one of 16 possible symbols (gray code mapping for 16-QPSK)
I = 2*symbols(:,1) - 1 + 2*(2*symbols(:,2) - 1);  % Mapping I component
Q = 2*symbols(:,3) - 1 + 2*(2*symbols(:,4) - 1);  % Mapping Q component

%% Modulation (RF Modulation)
fc = 10e3;               % Carrier frequency (Hz)
I_repeated = repelem(I, Ns);  % Repeat I for each sample
Q_repeated = repelem(Q, Ns);  % Repeat Q for each sample

t = (0:length(I_repeated)-1) * Ts; % Time vector for the entire signal
carrier_I = cos(2*pi*fc*t);  % I component carrier (cos)
carrier_Q = sin(2*pi*fc*t);  % Q component carrier (sin)

% Repeat I and Q components to match the oversampling rate


% Modulate I and Q components
tx_signal = I_repeated .* carrier_I - Q_repeated .* carrier_Q;  % RF modulated signal

%% AWGN Channel
for snr_idx = 1:length(Eb_N0_dB)
    Eb_N0 = 10^(Eb_N0_dB(snr_idx)/10);  % Convert dB to linear scale
    noise_power = (1/(2*Eb_N0)); % Noise power (baseband)
    noise = sqrt(noise_power) * randn(size(tx_signal)); % AWGN

    % Received signal with noise
    rx_signal = tx_signal + noise;

    %% Receiver Low-Pass Filtering (Matched Filter)
    h = ones(1, Ns);  % Simple rectangular matched filter
    rx_filtered = conv(rx_signal, h, 'same') / Ns; % Low-pass filtering and normalization

    %% Symbol Sampling & Decision
    sampled_symbols_I = rx_filtered(Ns:Ns:end);  % Sample at symbol centers for I component
    sampled_symbols_Q = rx_filtered(Ns/2:Ns:end); % Sample at symbol centers for Q component

    % Decision rule for demodulation
    demod_bits_I = sampled_symbols_I > 0;  % I channel: 1 if positive, 0 otherwise
    demod_bits_Q = sampled_symbols_Q > 0;  % Q channel: 1 if positive, 0 otherwise

    % Combine I and Q bits to get final demodulated bits
    demod_bits = zeros(1, N);
    demod_bits(1:4:end) = demod_bits_I;  % Assign demodulated I bits
    demod_bits(2:4:end) = demod_bits_Q;  % Assign demodulated Q bits

    %% BER Calculation
    CNT(snr_idx) = sum(bits ~= demod_bits);  % Count bit errors
    BER(snr_idx) = CNT(snr_idx) / N;       % Calculate Bit Error Rate
end

%% Plot BER vs Eb/N0
figure;
semilogy(Eb_N0_dB, BER, 'bo-', 'LineWidth', 2);
grid on;
xlabel('E_b/N_0 (dB)');
ylabel('Bit Error Rate (BER)');
title('16-QPSK BER Performance over AWGN (Realistic RF Modulation and Demodulation)');
legend('Simulated BER');
