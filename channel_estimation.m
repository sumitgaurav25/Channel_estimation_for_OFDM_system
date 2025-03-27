% Complete OFDM system with LS and MMSE channel estimation

% System parameters
N = 64;                 % Number of subcarriers
Ncp = 16;               % Cyclic prefix length
Np = 8;                 % Number of pilots
N_data = N - Np;        % Number of data subcarriers
mod_order = 16;         % Modulation order (16-QAM)
N_symbols = 100;        % Number of OFDM symbols
SNR_dB = 0:5:30;        % SNR range in dB

% Pilot arrangement (comb-type)
pilot_loc = 1:N/Np:N;   % Equally spaced pilot locations
data_loc = setdiff(1:N, pilot_loc); % Data subcarrier locations

% Pilot symbols (BPSK for simplicity)
pilots = (2*randi([0 1], 1, Np) - 1) + 1i*(2*randi([0 1], 1, Np) - 1);
pilots = pilots/sqrt(2); % Normalize power

% Generate random data symbols
data = randi([0 mod_order-1], N_data, N_symbols);

% Modulate data using QAM
mod_data = qammod(data, mod_order, 'UnitAveragePower', true);

% Initialize OFDM frame
ofdm_symbols = zeros(N, N_symbols);

% Insert pilots and data
for k = 1:N_symbols
    ofdm_symbols(pilot_loc, k) = pilots;
    ofdm_symbols(data_loc, k) = mod_data(:, k);
end

% Convert to time domain with IFFT
tx_signal = ifft(ofdm_symbols, N);

% Add cyclic prefix
tx_signal_cp = [tx_signal(end-Ncp+1:end, :); tx_signal];

% Channel parameters
max_delay = 10;         % Maximum delay spread (samples)
channel_taps = 4;       % Number of channel taps
h = (randn(1, channel_taps) + 1i*randn(1, channel_taps))/sqrt(2);
h = h/norm(h);          % Normalize channel power

% Frequency response of the channel
H_freq = fft(h, N);

% Initialize received signal
rx_signal_cp = zeros(size(tx_signal_cp));

% Channel convolution (time domain)
for k = 1:N_symbols
    rx_signal_cp(:, k) = conv(tx_signal_cp(:, k), h.', 'same');
end

% Initialize BER results
ber_ls = zeros(size(SNR_dB));
ber_mmse = zeros(size(SNR_dB));

% Channel autocorrelation matrix (corrected implementation)
R_hh = zeros(N, N);
for n = 1:N
    for m = 1:N
        R_hh(n,m) = sum(h.*conj(h).*exp(-1i*2*pi*(n-m)*(0:length(h)-1)/N));
    end
end

for snr_idx = 1:length(SNR_dB)
    % Convert SNR from dB to linear
    SNR_linear = 10^(SNR_dB(snr_idx)/10);
    
    % Calculate noise power
    signal_power = mean(abs(rx_signal_cp(:)).^2);
    noise_power = signal_power / SNR_linear;
    
    % Generate complex AWGN noise
    noise = sqrt(noise_power/2) * (randn(size(rx_signal_cp)) + 1i*randn(size(rx_signal_cp)));
    
    % Add noise to received signal
    rx_signal_noisy = rx_signal_cp + noise;
    
    % Remove cyclic prefix
    rx_signal = rx_signal_noisy(Ncp+1:end, :);
    
    % Convert to frequency domain with FFT
    rx_symbols = fft(rx_signal, N);
    
    % Extract received pilots
    rx_pilots = rx_symbols(pilot_loc, :);
    
    % LS Channel Estimation
    H_ls_pilots = rx_pilots ./ repmat(pilots.', 1, N_symbols);
    
    % Interpolate to estimate channel at all subcarriers
    H_ls = zeros(N, N_symbols);
    for k = 1:N_symbols
        H_ls(:, k) = interp1(pilot_loc, H_ls_pilots(:, k), 1:N, 'spline');
    end
    
    % MMSE Channel Estimation
    SNR_linear_pilots = SNR_linear * N/Np; % Pilot SNR is boosted
    
    % MMSE estimation matrix
    W_mmse = R_hh / (R_hh + (1/SNR_linear_pilots)*eye(N));
    
    % Apply MMSE estimation
    H_mmse = zeros(N, N_symbols);
    for k = 1:N_symbols
        H_mmse(:, k) = W_mmse * H_ls(:, k);
    end
    
    % Equalization and data detection
    rx_data_ls = zeros(N_data, N_symbols);
    rx_data_mmse = zeros(N_data, N_symbols);
    
    for k = 1:N_symbols
        % LS equalization
        rx_data_ls(:, k) = rx_symbols(data_loc, k) ./ H_ls(data_loc, k);
        
        % MMSE equalization
        rx_data_mmse(:, k) = rx_symbols(data_loc, k) ./ H_mmse(data_loc, k);
    end
    
    % Demodulate data
    demod_data_ls = qamdemod(rx_data_ls, mod_order, 'UnitAveragePower', true);
    demod_data_mmse = qamdemod(rx_data_mmse, mod_order, 'UnitAveragePower', true);
    
    % Calculate BER
    [~, ber_ls(snr_idx)] = biterr(data(:), demod_data_ls(:));
    [~, ber_mmse(snr_idx)] = biterr(data(:), demod_data_mmse(:));
end

% Plot actual vs estimated channel (for one symbol)
figure;
subplot(2,1,1);
plot(1:N, abs(H_freq), 'b', 1:N, abs(H_ls(:,1)), 'r--');
title('LS Channel Estimation');
legend('Actual Channel', 'Estimated Channel');
xlabel('Subcarrier Index'); ylabel('Magnitude');

subplot(2,1,2);
plot(1:N, abs(H_freq), 'b', 1:N, abs(H_mmse(:,1)), 'g--');
title('MMSE Channel Estimation');
legend('Actual Channel', 'Estimated Channel');
xlabel('Subcarrier Index'); ylabel('Magnitude');

% Plot BER vs SNR
figure;
semilogy(SNR_dB, ber_ls, 'r-o', SNR_dB, ber_mmse, 'b-s', 'LineWidth', 2);
grid on;
xlabel('SNR (dB)'); ylabel('Bit Error Rate (BER)');
title('BER Performance of LS vs MMSE Channel Estimation');
legend('LS Estimation', 'MMSE Estimation');