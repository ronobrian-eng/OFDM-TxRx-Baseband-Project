% ========================================================================
% Baseband OFDM Transmitter / Receiver (Custom 16-QAM, AWGN, 3-Tap Channel)
% Author: Brian Rono
%
% - 16-QAM on N subcarriers
% - IFFT + cyclic prefix at TX
% - Static 3-tap multipath channel + AWGN (manual)
% - CP removal + FFT + one-tap equaliser at RX
% - BER vs SNR curve
%
% Figures:
%   1) Magnitude of one OFDM symbol (with CP)
%   2) 16-QAM constellation: TX vs equalised RX (mid SNR)
%   3) Channel frequency response |H[k]|
%   4) BER vs SNR (semilogy)
%
% Uses only base MATLAB (no toolboxes required).
% ========================================================================

clear; close all; clc;

% ------------------------ 1. System Parameters --------------------------
Nsc   = 64;        % number of subcarriers
Ncp   = 16;        % cyclic prefix length (samples)
M     = 16;        % 16-QAM
k_mod = 4;         % bits per symbol for 16-QAM

Nsym  = 1000;      % OFDM symbols used for BER statistics

SNRdB_vec = 0:5:30;            % SNR points for BER curve (dB)
Nsamples_per_ofdm = Nsc + Ncp; % samples per OFDM symbol including CP

fprintf('OFDM: %d subcarriers, CP = %d, %d-QAM, %d symbols\n', ...
    Nsc, Ncp, M, Nsym);

% ------------------------ 2. Simple Multipath Channel -------------------
% Static 3-tap baseband channel (can be adjusted)
h = [0.9; 0.4*exp(1j*pi/4); 0.2*exp(-1j*pi/3)];
Lh = length(h);

if Ncp < (Lh-1)
    warning('CP is shorter than channel length; ISI will occur.');
end

% Frequency response per subcarrier (for equalisation and plotting)
Hk = fft(h, Nsc).';     % Nsc x 1

% ------------------------ 3. Generate Random 16-QAM Symbols ------------
% Total QAM symbols: Nsc subcarriers * Nsym OFDM symbols
N_qam  = Nsc * Nsym;
N_bits = N_qam * k_mod;

% Random bit stream
tx_bits_vec = randi([0 1], N_bits, 1);

% Group bits into 4 per symbol
tx_bits = reshape(tx_bits_vec, k_mod, N_qam).';  % N_qam x 4

% Map bits to 16-QAM points (custom mapper)
data_qam = bits_to_16qam(tx_bits);               % N_qam x 1 (complex)

% Arrange as Nsc x Nsym in frequency domain
X = reshape(data_qam, Nsc, Nsym);                % Nsc x Nsym

% ------------------------ 4. OFDM Modulator (IFFT + CP) -----------------
% Time-domain OFDM symbols (without CP)
x_noCP = ifft(X, Nsc, 1);                        % Nsc x Nsym

% Add cyclic prefix
x_cp = [x_noCP(end-Ncp+1:end,:); x_noCP];        % (Nsc+Ncp) x Nsym

% Serialised transmit signal (all OFDM symbols)
tx_signal = x_cp(:); %#ok<NASGU>

% ------------------------ 5. Plot One OFDM Symbol (Time Domain) ---------
figure;
example_sym_idx = 1;
example_sym = x_cp(:, example_sym_idx);
stem(0:Nsamples_per_ofdm-1, abs(example_sym), 'filled');
xlabel('Sample index');
ylabel('|x[n]|');
title('Magnitude of One OFDM Symbol (with Cyclic Prefix)');
grid on;

% ------------------------ 6. Channel + Noise + Detection ----------------
BER = zeros(size(SNRdB_vec));

% For constellation snapshot at a particular SNR
SNRdB_snapshot = 15;
constellation_tx_snapshot = [];
constellation_rx_snapshot = [];

for ii = 1:numel(SNRdB_vec)
    SNRdB = SNRdB_vec(ii);

    % ----- Channel: apply convolution per OFDM symbol --------------------
    rx_cp = zeros(Nsamples_per_ofdm, Nsym);
    for n = 1:Nsym
        this_sym = x_cp(:, n);            % (Nsc+Ncp) x 1
        y_conv   = conv(this_sym, h);     % length Nsamples_per_ofdm + Lh - 1
        rx_cp(:, n) = y_conv(1:Nsamples_per_ofdm);  % truncate to symbol length
    end

    % Serialise the channel output
    rx_serial = rx_cp(:);

    % ----- Add AWGN (manual) ---------------------------------------------
    % Signal power
    Ps = mean(abs(rx_serial).^2);
    % Desired SNR in linear scale
    SNR_lin = 10^(SNRdB/10);
    % Noise power
    Pn = Ps / SNR_lin;
    % Complex Gaussian noise
    noise = sqrt(Pn/2) * (randn(size(rx_serial)) + 1j*randn(size(rx_serial)));

    rx_noisy = rx_serial + noise;

    % Reshape back to (Nsc+Ncp) x Nsym
    rx_cp_noisy = reshape(rx_noisy, Nsamples_per_ofdm, Nsym);

    % ----- Remove CP and perform FFT ------------------------------------
    y_noCP = rx_cp_noisy(Ncp+1:end, :);      % Nsc x Nsym
    Y = fft(y_noCP, Nsc, 1);                 % frequency-domain with channel

    % ----- One-tap equalisation on each subcarrier ----------------------
    Y_eq = Y ./ repmat(Hk, 1, Nsym);        % Nsc x Nsym

    % Collect equalised QAM symbols
    data_qam_hat = Y_eq(:);                 % N_qam x 1

    % ----- 16-QAM Demapping ---------------------------------------------
    rx_bits = qam16_to_bits(data_qam_hat);  % N_qam x 4

    rx_bits_vec = rx_bits.';                % 4 x N_qam
    rx_bits_vec = rx_bits_vec(:);           % (N_qam*4) x 1

    % ----- Bit errors ---------------------------------------------------
    numErr  = sum(tx_bits_vec ~= rx_bits_vec);
    numBits = numel(tx_bits_vec);
    BER(ii) = numErr / numBits;

    % Save constellation snapshot at chosen SNR
    if abs(SNRdB - SNRdB_snapshot) < 1e-9
        constellation_tx_snapshot = X(:);          % original 16-QAM symbols
        constellation_rx_snapshot = data_qam_hat;  % after equalisation
    end

    fprintf('SNR = %2d dB -> BER = %.3e\n', SNRdB, BER(ii));
end

% ------------------------ 7. Constellations (TX vs RX) ------------------
if ~isempty(constellation_tx_snapshot)
    max_points_plot = 3000;
    tx_plot = constellation_tx_snapshot(1:min(end,max_points_plot));
    rx_plot = constellation_rx_snapshot(1:min(end,max_points_plot));

    figure;
    subplot(1,2,1);
    plot(real(tx_plot), imag(tx_plot), '.', 'MarkerSize', 6);
    axis equal; grid on;
    xlabel('In-Phase'); ylabel('Quadrature');
    title('Transmitted 16-QAM Constellation');

    subplot(1,2,2);
    plot(real(rx_plot), imag(rx_plot), '.', 'MarkerSize', 6);
    axis equal; grid on;
    xlabel('In-Phase'); ylabel('Quadrature');
    title(sprintf('Received (Equalised) 16-QAM, SNR = %d dB', SNRdB_snapshot));
end

% ------------------------ 8. Channel Frequency Response -----------------
figure;
stem(0:Nsc-1, abs(Hk), 'filled');
xlabel('Subcarrier index k');
ylabel('|H[k]|');
title('Magnitude of Channel Frequency Response |H[k]|');
grid on;

% ------------------------ 9. BER vs SNR Curve ---------------------------
figure;
semilogy(SNRdB_vec, BER, '-o', 'LineWidth', 1.5);
grid on;
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs SNR for 16-QAM OFDM over 3-Tap Channel');
ylim([1e-5 1]);
xlim([min(SNRdB_vec) max(SNRdB_vec)]);

disp('OFDM simulation finished.');

% ========================================================================
% Helper functions
% ========================================================================

function s = bits_to_16qam(bits)
%BITS_TO_16QAM Map 4-bit rows to 16-QAM constellation points.
% bits: N x 4 (each row is [b1 b2 b3 b4])

    bits = logical(bits);
    N = size(bits,1);

    % Map first 2 bits -> I, last 2 bits -> Q
    bI = bits(:,1:2);   % N x 2
    bQ = bits(:,3:4);   % N x 2

    % Natural mapping: 00 -> -3, 01 -> -1, 10 -> +1, 11 -> +3
    lut = [-3; -1; 1; 3];

    idxI = bI(:,1)*2 + bI(:,2);       % 0..3
    idxQ = bQ(:,1)*2 + bQ(:,2);       % 0..3

    I = lut(idxI+1);
    Q = lut(idxQ+1);

    % Normalise so average symbol power is 1
    % For 16-QAM with levels {-3,-1,1,3}, Es_avg = 10
    s = (I + 1j*Q) / sqrt(10);

end

function bits = qam16_to_bits(s)
%QAM16_TO_BITS Demap 16-QAM symbols back to 4-bit rows.
% s: N x 1 complex symbols

    % Undo normalisation
    s = s * sqrt(10);

    I = real(s);
    Q = imag(s);

    % Quantise I and Q to nearest in {-3,-1,1,3}
    levels = [-3 -1 1 3];

    Iq = quantise_to_levels(I, levels);
    Qq = quantise_to_levels(Q, levels);

    % Convert quantised amplitudes back to 2-bit indices
    [~, idxI] = ismember(Iq, levels);
    [~, idxQ] = ismember(Qq, levels);

    idxI = idxI - 1;    % 0..3
    idxQ = idxQ - 1;    % 0..3

    % Convert indices to bits: 0..3 -> 2-bit binary
    N = numel(s);
    bitsI = zeros(N, 2);
    bitsQ = zeros(N, 2);

    bitsI(:,1) = floor(idxI/2);
    bitsI(:,2) = mod(idxI,2);

    bitsQ(:,1) = floor(idxQ/2);
    bitsQ(:,2) = mod(idxQ,2);

    bits = [bitsI bitsQ];    % N x 4
end

function xq = quantise_to_levels(x, levels)
%QUANTISE_TO_LEVELS Quantise real vector x to nearest value in "levels".
    xq = zeros(size(x));
    for k = 1:numel(x)
        [~, idx] = min(abs(x(k) - levels));
        xq(k) = levels(idx);
    end
end
