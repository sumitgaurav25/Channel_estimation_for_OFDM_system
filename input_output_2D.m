clear all; close all; clc;

%% Load and Preprocess Image
img = imread('lena.bmp'); % Load image
img_gray = rgb2gray(img); % Convert to grayscale
img_down = imresize(img_gray, 0.25); % Downsample by factor of 4

%% Convert Image to Binary Stream
[rows, cols] = size(img_down);
bits_per_pixel = 8;
image_bits = reshape(de2bi(img_down, bits_per_pixel)', 1, []); % Convert pixels to binary

%% OFDM System Parameters
N = 64; Ncp = 16; mod_order = 16;
bits_per_symbol = log2(mod_order);
num_symbols = ceil(length(image_bits) / (bits_per_symbol * N));

%% Zero Padding
padded_bits = [image_bits, zeros(1, N * num_symbols * bits_per_symbol - length(image_bits))];

%% QAM Modulation
reshaped_bits = reshape(padded_bits, bits_per_symbol, [])';
qam_symbols = qammod(bi2de(reshaped_bits), mod_order, 'UnitAveragePower', true);
ofdm_symbols = reshape(qam_symbols, N, num_symbols);

%% IFFT and Cyclic Prefix
tx_time = ifft(ofdm_symbols, N);
tx_with_cp = [tx_time(end-Ncp+1:end, :); tx_time];
tx_signal = tx_with_cp(:).';

%% Channel Model and Noise
channel_taps = 4;
h = (randn(1, channel_taps) + 1i*randn(1, channel_taps)) / sqrt(2);
h = h / norm(h);
rx_signal = conv(tx_signal, h, 'same');
SNR_dB = 20;
noise_power = mean(abs(rx_signal).^2) / (10^(SNR_dB/10));
noise = sqrt(noise_power/2) * (randn(size(rx_signal)) + 1i*randn(size(rx_signal)));
rx_signal_noisy = rx_signal + noise;

%% Receiver Processing
rx_reshaped = reshape(rx_signal_noisy, N+Ncp, num_symbols);
rx_no_cp = rx_reshaped(Ncp+1:end, :);
rx_freq = fft(rx_no_cp, N);

%% LS Channel Estimation
H_ls = mean(rx_freq ./ ofdm_symbols, 2);
%rx_data_ls = rx_freq ./ H_ls;
%rx_serial_ls = rx_data_ls(:).';

rx_data_ls = rx_freq ./ H_ls; 
rx_serial_ls = rx_data_ls(:).';


% demod_symbols_ls = qamdemod(rx_serial_ls, mod_order, 'UnitAveragePower', true);
% output_bits_ls = reshape(de2bi(demod_symbols_ls, bits_per_symbol)', 1, []);
% output_bits_ls = output_bits_ls(1:length(image_bits));

demod_symbols_ls = qamdemod(rx_serial_ls, mod_order, 'UnitAveragePower', true);
output_bits_ls = reshape(de2bi(demod_symbols_ls, bits_per_symbol)', 1, []);
output_bits_ls = output_bits_ls(1:length(image_bits)); % Trim to match original length


%% MMSE Channel Estimation
SNR_linear = 10^(SNR_dB/10);
H_mmse = H_ls ./ (1 + (1/SNR_linear));
rx_data_mmse = rx_freq ./ H_mmse;
rx_serial_mmse = rx_data_mmse(:).';

demod_symbols_mmse = qamdemod(rx_serial_mmse, mod_order, 'UnitAveragePower', true);
output_bits_mmse = reshape(de2bi(demod_symbols_mmse, bits_per_symbol)', 1, []);
output_bits_mmse = output_bits_mmse(1:length(image_bits));

%% Reconstruct Images
received_img_ls = reshape(uint8(bi2de(reshape(output_bits_ls, bits_per_pixel, [])')), rows, cols);
received_img_mmse = reshape(uint8(bi2de(reshape(output_bits_mmse, bits_per_pixel, [])')), rows, cols);

%% Display Images
figure;
subplot(1,3,1); imshow(img_down); title('Original Image');
subplot(1,3,2); imshow(received_img_ls); title('Received Image (LS)');
subplot(1,3,3); imshow(received_img_mmse); title('Received Image (MMSE)');

%% Plot Constellation Diagrams
figure;
subplot(1,2,1); plot(rx_serial_ls, 'r.'); title('Constellation (LS)');
subplot(1,2,2); plot(rx_serial_mmse, 'b.'); title('Constellation (MMSE)');

%% Bit Error Rate (BER) Calculation
ber_ls = sum(image_bits ~= output_bits_ls) / length(image_bits);
ber_mmse = sum(image_bits ~= output_bits_mmse) / length(image_bits);
fprintf('BER with LS Estimation: %.4f\n', ber_ls);
fprintf('BER with MMSE Estimation: %.4f\n', ber_mmse);

%% Bit Error Comparison

figure;
stem(abs(double(image_bits(1:100)) - double(output_bits_ls(1:100))), 'r', 'filled');
hold on;
stem(abs(double(image_bits(1:100)) - double(output_bits_mmse(1:100))), 'g', 'filled');
%hold off;
title('Bit Error Comparison'); xlabel('Bit Position'); ylabel('Error (1=wrong)');
legend('LS Errors', 'MMSE Errors');


subplot(2,1,1);
%stem(abs(double(image_bits(1:compare_bits)) - double(output_bits_ls(1:compare_bits))), 'r', 'filled');

stem(abs(double(image_bits(1:100)) - double(output_bits_ls(1:100))), 'r', 'filled');
title('LS Bit Errors');
ylabel('Error (1=wrong)');
xlabel('Bit Position');
%ylim([0 1.2]);

subplot(2,1,2);
%stem(abs(double(image_bits(1:compare_bits)) - double(output_bits_mmse(1:compare_bits))), 'b', 'filled');

stem(abs(double(image_bits(1:100)) - double(output_bits_mmse(1:100))), 'g', 'filled');
title('MMSE Bit Errors');
ylabel('Error (1=wrong)');
xlabel('Bit Position');
%ylim([0 1.2]);