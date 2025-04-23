clc; clear;

%% Parameters
ovs = 10;
N = 1e6;         % Number of bits
fc = 1e3;        % Carrier frequency (Hz)
fs = ovs*fc;      % Sampling frequency (Hz) (should be >> fc)
Ts = 1/fs;       % Sampling period
Eb_N0_dB = -2:0.2:2; % SNR values in dB
T_symbol = 1/fc; % Symbol duration (1 cycle per symbol)
Ns = fs / fc;    % Samples per symbol

%% Generate Random Bits
bits = randi([0 1], 1, N);  % Random 0s and 1s

%% BPSK Modulation (RF Domain)
t = 0:Ts:T_symbol-Ts;       % Time vector for one symbol
carrier_cos = cos(2*pi*fc*t);  % Carrier (I)
carrier_sin = sin(2*pi*fc*t);  % Carrier (Q) (unused in BPSK)

tx_signal = zeros(1,ovs*N);  % Initialize transmitted signal
for i = 1:N
    if bits(i) == 0
        mod_symbol = -carrier_cos;  % BPSK maps 0 -> -cos
    else
        mod_symbol = carrier_cos;   % BPSK maps 1 -> +cos
    end
    tx_signal((i-1)*ovs+1:i*ovs) =  mod_symbol;  % Append modulated symbol
end

%% AWGN Channel
BER = zeros(1, length(Eb_N0_dB)); % Store Bit Error Rate
errors = BER;
for snr_idx = 1:length(Eb_N0_dB)
    Eb_N0 = 10^(Eb_N0_dB(snr_idx)/10);  % Convert dB to linear scale
    noise_power = (1/(2*Eb_N0)); % Noise power (assuming unit energy)
    noise = sqrt(noise_power) * randn(size(tx_signal)); % AWGN

    rx_signal = tx_signal + noise;  % Received signal with noise

    %% Coherent Demodulation
    demod_bits = zeros(1, N);  
    for i = 1:N
        % Extract the symbol
        rx_symbol = rx_signal((i-1)*Ns + (1:Ns)); 

        % Multiply by cosine and integrate (dot product)
        I_component = sum(rx_symbol .* carrier_cos);
%     plot([rx_symbol; carrier_cos]')

        % Decision rule: If positive, bit = 1; else, bit = 0
        if I_component > 0
            demod_bits(i) = 1;
        else
            demod_bits(i) =  0;
        end
    end

    %% BER Calculation
    errors(snr_idx) = sum(bits ~= demod_bits);
    BER(snr_idx) = errors(snr_idx) / N;
end

%% Plot BER vs Eb/N0
figure;
semilogy(Eb_N0_dB, BER, 'ro-', 'LineWidth', 2);
grid on;
xlabel('E_b/N_0 (dB)');
ylabel('Bit Error Rate (BER)');
title('BPSK BER Performance over AWGN');
