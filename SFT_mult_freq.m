%% SFT Project:
% Implement the SFT algorithm in MATLAB. Make comparisons with the built-in 
% FFT algorithm of MATLAB in terms of running time as the length of data 
% increases, for different sparsity cases


%% create signal
% need to be able to make signals with N = length,k = number of peaks, 
% p = spacing of peaks, epsilon[j] = noise function

K = 1; % amount of peaks
A = 1; % amplitude of peaks
% use the folloing 3 lines for gaussian-like magnitudes
%Amax = K; % peak magnitude for central sample
%x = 1:K;
%A = K/Amax + (Amax - K/Amax) * (1 - abs(x - (K+1)/2)/((K+1)/2-1));
k0 = 35; % center frequency index, 1-index
N = 70; % signal length. MUST BE EVEN
delta_w = 0.7; % peak spacing, (0 to 1]. approach 0 for closer spacing
sampling_rate = 0.25; %amount of samples to use for fourier coeff estimation, (0,1]. Higher value is more accuracy but slower

n  = 0:N-1; % time indices


k_list = (1:N/K:N) + ((k0-512) + N/(2*K));
half_k_range = round(((k_list(end) - k_list(1)) * delta_w) / 2);   % half of total range
k_list = linspace(k0 - half_k_range, k0 + half_k_range, K);
k_list = round(k_list);
k = zeros(1,N);
k(k_list) = A;
fprintf("setting peaks at %d\n", k_list)

%plot freq peaks
indices = 1:length(k);
nonzero_idx = k ~= 0;
figure;
stem(indices(nonzero_idx), k(nonzero_idx), 'o');
% use following two lines for x ticks for only peaks
%ax = gca; 
%ax.XTick = indices(nonzero_idx);
xlim([0 N]);
xlabel('Index');
ylabel('Magnitude');
title('Sparse Frequency Peaks');
grid off;

% Create time signal
S = exp(2*pi*1i * (k_list(:) * n) / N);
f_n = sum(S, 1);    % 1-by-N signal
f_orig = f_n;
f_k_hat = zeros(1,N); %assume we are taking an N-pt DFT of f_n

%sampling for fourier coeff estimation
L = floor(sampling_rate * N);

%% Create filter
a = 0.5; %transition band width
B = K*4; %number of freq bins
m_pass = floor((1-a)*N / (2*B)); %passband width 
m_stop = floor(N/(2*B));
k_passband = m_pass+2:m_stop+1; %freq index variable, same length as signal for an N-pt SFT
gamma = 1; %smoothness of gaussian taper 
G_k_passband =  exp(-1 * (abs(k_passband) - m_pass).^2 ./ (2*gamma^2));
G_k = zeros(1,N);
G_k(m_pass+2:m_stop+1) = G_k_passband;
G_k(1:m_pass+1) = 1;
g_n = ifft(G_k);

%plot freq response (gaussian)
figure;
stem(abs(G_k))
xlabel('k');
ylabel('Magnitude');
title('Frequency Filter');
grid off;

%plot time response
figure;
stem(abs(g_n))
xlabel('n');
ylabel('Magnitude');
title('Time Filter');
grid off;

%downsample filter appropriately
g_prime_n = g_n(0:floor(N/B):B-1);

%plot downsampled time response
figure;
stem(abs(g_prime_n))
xlabel('n');
ylabel('Magnitude');
title('Downsampled Time Filter');
grid off;

%% Main loop: randomly permute, find recovered frequencies, then peel
while true
    %% Permute signal
    n = 0:N-1;
    sigma = randi([1, floor((N-1)/2)]) * 2 - 1; %sigma can be any odd value [1,N-1]
    a = randi([0, N-1]);
    b = randi([0, N-1]);
    idx = mod(sigma*n + b, N) + 1;
    idx_no_tshift = mod(sigma*n, N) + 1;
    f_sigma = f_n(idx);
    f_notshift = f_n(idx_no_tshift);
    f_perm = f_sigma .* exp(-1j*2*pi*a*n/N);

    %calculate sigma inverse mod N for unpermute
    [~, x_temp, ~] = gcd(sigma, N);
    sigma_inv = mod(x_temp,N);

    %% Implement SFT I. Identify Frequency Peaks

    % Use random binning with modulo aliasing and CRT

    % Random binning

    % Aliasing and CRT
    % calculate padded signal of length with specific coprime factors.
    % subsampling rates MUST be > k but << N
    subsampling_lengths = [2,5,7]; 
    subsampling_rates = [N/subsampling_lengths(1), N/subsampling_lengths(2), N/subsampling_lengths(3)];

    % find subsampled DFTs
    f_1 = fft(f_perm(1:subsampling_rates(1):end));
    f_2 = fft(f_perm(1:subsampling_rates(2):end));
    f_3 = fft(f_perm(1:subsampling_rates(3):end));

    % debug
    [B_1, ~] = sort(abs(f_1), 'descend');
    [B_2, ~] = sort(abs(f_2), 'descend');
    [B_3, ~] = sort(abs(f_3), 'descend');

    %threshold for noise robustness
    thresh1 = 0.85 * max(abs(f_1));
    thresh2 = 0.85 * max(abs(f_2));
    thresh3 = 0.85 * max(abs(f_3));

    %find peaks
    peaks1 = find(abs(f_1) > thresh1);
    peaks2 = find(abs(f_2) > thresh2);
    peaks3 = find(abs(f_3) > thresh3);

    if isempty(peaks1) || isempty(peaks2) || isempty(peaks3)
        disp("no peaks found")
        break
    end

    k_perm = crt([peaks1(1)-1, peaks2(1)-1, peaks3(1)-1],subsampling_lengths); 
    %convert from 1- to 0- index for CRT above, then back to 1-index

    k_hat = mod((k_perm+a)*sigma_inv, N);

    %fourier coefficient estimation
    samp_idxs = randi([1 N], 1, L);
    f_samp = f_n(samp_idxs);
    A_hat = 1/L * sum(f_samp.*exp(-2j*pi*k_hat*samp_idxs/N));
    if max(f_k_hat) ~= 0 && abs(A_hat) < 0.1 * abs(max(f_k_hat))
        disp("peaks below threshold magnitude")
        break
    end

    f_k_hat(k_hat) = A_hat;
    fprintf('recovered peaks at %d with magnitude %d\n', k_hat, abs(A_hat))

    %peel
    disp('peeling...')

    %debug:     peel_contrib = max(fft(A_hat .* exp(-2j*pi*n*k_hat/N)));
    f_n = f_n + A_hat .* exp(-2j*pi*n*k_hat/N);

    tempval = 'debug';
end

disp('done.')
%disp(f_k_hat)