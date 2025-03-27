clear all; close all; clc;

%% System Parameters
N = 64;                 % Number of subcarriers
Ncp = 16;               % Cyclic prefix length
mod_order = 16;         % Modulation order (16-QAM)
num_symbols = 10;       % Number of OFDM symbols to transmit
SNR_dB = 20;           % SNR in dB

%% 1. Generate Random Binary Input Sequence
bits_per_symbol = log2(mod_order); % Bits per QAM symbol
total_bits = bits_per_symbol * N * num_symbols; % Total bits needed
input_bits = randi([0 1], 1, total_bits); % Random binary sequence

%% 2. QAM Modulation
% Reshape bits into groups for QAM symbols (each column has bits_per_symbol bits)
reshaped_bits = reshape(input_bits, bits_per_symbol, [])';
qam_symbols = qammod(bi2de(reshaped_bits), mod_order, 'UnitAveragePower', true);

%% 3. Serial-to-Parallel Conversion
% Reshape into OFDM symbols (each column is one OFDM symbol)
ofdm_symbols = reshape(qam_symbols, N, num_symbols);

%% 4. IFFT (Time Domain Conversion)
tx_time = ifft(ofdm_symbols, N);

%% 5. Parallel-to-Serial Conversion
tx_serial = tx_time(:).'; % Convert to serial stream

%% 6. Add Cyclic Prefix
% Reshape to add CP to each symbol
tx_with_cp = reshape(tx_serial, N, num_symbols);
tx_with_cp = [tx_with_cp(end-Ncp+1:end, :); tx_with_cp];
tx_signal = tx_with_cp(:).'; % Final transmitted signal

%% Channel Model
% Create frequency-selective channel
channel_taps = 4; % Number of channel taps
h = (randn(1, channel_taps) + 1i*randn(1, channel_taps))/sqrt(2);
h = h/norm(h); % Normalize channel power

% Channel convolution
rx_signal = conv(tx_signal, h, 'same');

% Add AWGN noise
signal_power = mean(abs(rx_signal).^2);
noise_power = signal_power / (10^(SNR_dB/10));
noise = sqrt(noise_power/2)*(randn(size(rx_signal)) + 1i*randn(size(rx_signal)));
rx_signal_noisy = rx_signal + noise;

%% Receiver Processing (Common steps before channel estimation)
% 1. Remove Cyclic Prefix
rx_reshaped = reshape(rx_signal_noisy, N+Ncp, num_symbols);
rx_no_cp = rx_reshaped(Ncp+1:end, :);

% 2. Serial-to-Parallel Conversion
rx_parallel = rx_no_cp;

% 3. FFT (Frequency Domain Conversion)
rx_freq = fft(rx_parallel, N);

%% Channel Estimation (Using Pilots)
% Pilot parameters
Np = 8; % Number of pilots
pilot_loc = 1:N/Np:N; % Equally spaced pilot locations
data_loc = setdiff(1:N, pilot_loc); % Data subcarrier locations

% Known pilot symbols (BPSK)
pilots = (2*randi([0 1], 1, Np) - 1) + 1i*(2*randi([0 1], 1, Np) - 1);
pilots = pilots/sqrt(2); % Normalize power

% Insert pilots into transmitted signal (for estimation)
tx_pilots = zeros(N, num_symbols);
tx_pilots(pilot_loc, :) = repmat(pilots.', 1, num_symbols);

%% LS Channel Estimation
H_ls_pilots = rx_freq(pilot_loc, :) ./ tx_pilots(pilot_loc, :);

% Interpolate to estimate channel at all subcarriers
H_ls = zeros(N, num_symbols);
for k = 1:num_symbols
    H_ls(:, k) = interp1(pilot_loc, H_ls_pilots(:, k), 1:N, 'spline');
end

%% MMSE Channel Estimation
% Channel autocorrelation matrix
R_hh = zeros(N, N);
for n = 1:N
    for m = 1:N
        R_hh(n,m) = sum(h.*conj(h).*exp(-1i*2*pi*(n-m)*(0:length(h)-1)/N));
    end
end

SNR_linear = 10^(SNR_dB/10);
SNR_linear_pilots = SNR_linear * N/Np; % Pilot SNR is boosted

% MMSE estimation matrix
W_mmse = R_hh / (R_hh + (1/SNR_linear_pilots)*eye(N));

% Apply MMSE estimation
H_mmse = zeros(N, num_symbols);
for k = 1:num_symbols
    H_mmse(:, k) = W_mmse * H_ls(:, k);
end

%% Equalization and Data Detection
% LS Equalization
rx_data_ls = rx_freq(data_loc, :) ./ H_ls(data_loc, :);

% MMSE Equalization
rx_data_mmse = rx_freq(data_loc, :) ./ H_mmse(data_loc, :);

%% Parallel-to-Serial Conversion
rx_serial_ls = rx_data_ls(:).';
rx_serial_mmse = rx_data_mmse(:).';

%% QAM Demodulation
% LS Demodulation
demod_symbols_ls = qamdemod(rx_serial_ls, mod_order, 'UnitAveragePower', true);
output_bits_ls = reshape(de2bi(demod_symbols_ls, bits_per_symbol)', 1, []);

% MMSE Demodulation
demod_symbols_mmse = qamdemod(rx_serial_mmse, mod_order, 'UnitAveragePower', true);
output_bits_mmse = reshape(de2bi(demod_symbols_mmse, bits_per_symbol)', 1, []);

%% Results
disp('Input Bits (First 40):');
disp(num2str(input_bits(1:40)));

disp('Output Bits with LS Estimation (First 40):');
disp(num2str(output_bits_ls(1:40)));

disp('Output Bits with MMSE Estimation (First 40):');
disp(num2str(output_bits_mmse(1:40)));

% Calculate BER
%ber_ls = sum(input_bits ~= output_bits_ls)/length(input_bits);

min_len = min(length(input_bits), length(output_bits_ls));
ber_ls = sum(input_bits(1:min_len) ~= output_bits_ls(1:min_len)) / min_len;

min_len = min(length(input_bits), length(output_bits_mmse));
ber_mmse = sum(input_bits(1:min_len) ~= output_bits_mmse(1:min_len)) / min_len;

fprintf('\nBER with LS Estimation: %.4f\n', ber_ls);
fprintf('BER with MMSE Estimation: %.4f\n', ber_mmse);

%% Plotting
% Plot channel estimates
figure;
subplot(2,1,1);
plot(1:N, abs(fft(h, N)), 'b', 1:N, abs(H_ls(:,1)), 'r--');
title('Channel Estimation Comparison (LS)');
legend('Actual Channel', 'Estimated Channel');
xlabel('Subcarrier Index'); ylabel('Magnitude');
grid on;

subplot(2,1,2);
plot(1:N, abs(fft(h, N)), 'b', 1:N, abs(H_mmse(:,1)), 'g--');
title('Channel Estimation Comparison (MMSE)');
legend('Actual Channel', 'Estimated Channel');
xlabel('Subcarrier Index'); ylabel('Magnitude');
grid on;

% Plot constellation diagrams
figure;
subplot(1,2,1);
plot(rx_serial_ls, 'r.');
title('Received Constellation (LS Estimation)');
xlabel('In-Phase'); ylabel('Quadrature');
axis square; grid on;

subplot(1,2,2);
plot(rx_serial_mmse, 'b.');
title('Received Constellation (MMSE Estimation)');
xlabel('In-Phase'); ylabel('Quadrature');
axis square; grid on;

% Plot bit error comparison
figure;
stem(1:40, input_bits(1:40), 'b', 'LineWidth', 2, 'DisplayName', 'Input');
hold on;
stem(1:40, output_bits_ls(1:40), 'r', 'DisplayName', 'LS Output');
stem(1:40, output_bits_mmse(1:40), 'g', 'DisplayName', 'MMSE Output');
hold off;
title('Bit Error Comparison (First 40 Bits)');
xlabel('Bit Index'); ylabel('Bit Value');
legend;
grid on;

% Create error difference plot
figure;
subplot(2,1,1);
stem(abs(input_bits(1:40) - output_bits_ls(1:40)), 'r', 'filled');
title('LS Estimation Error Pattern');
ylabel('Error (1=wrong)');
xlabel('Bit Position');
ylim([0 1.2]);
grid on;

subplot(2,1,2);
stem(abs(input_bits(1:40) - output_bits_mmse(1:40)), 'g', 'filled');
title('MMSE Estimation Error Pattern');
ylabel('Error (1=wrong)');
xlabel('Bit Position');
ylim([0 1.2]);
grid on;