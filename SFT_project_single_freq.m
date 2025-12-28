%% SFT Project:
% Implement the SFT algorithm in MATLAB. Make comparisons with the built-in 
% FFT algorithm of MATLAB in terms of running time as the length of data 
% increases, for different sparsity cases


%% create signal
% need to be able to make signals with N = length,k = number of peaks, 
% p = spacing of peaks, epsilon[j] = noise function

K = 1; % amount of peaks
A = 1;
% use the folloing 3 lines for gaussian-like magnitudes
%Amax = K; % peak magnitude for central sample
%x = 1:K;
%A = K/Amax + (Amax - K/Amax) * (1 - abs(x - (K+1)/2)/((K+1)/2-1));
k0 = 59; % center frequency index
N = 70; % signal length
delta_w = 0.6; % peak spacing, (0 to 1]. approach 0 for closer spacing
n  = 0:N-1; % time indices

k_list = (1:N/K:N) + ((k0-512) + N/(2*K));
half_k_range = round(((k_list(end) - k_list(1)) * delta_w) / 2);   % half of total range
k_list = linspace(k0 - half_k_range, k0 + half_k_range, K);
k_list = round(k_list);
k = zeros(1,N);
k(k_list) = A;

%plot freq peaks
indices = 1:length(k);
nonzero_idx = k ~= 0;
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
f = sum(S, 1);    % 1-by-N signal

%% Implement SFT I. Identify Frequency Peaks

% Use random binning with modulo aliasing and CRT

% Random binning

% Aliasing and CRT
% calculate padded signal of length with specific coprime factors.
% subsampling rates MUST be > k but << N
subsampling_lengths = [2,5,7]; 
subsampling_rates = [N/subsampling_lengths(1), N/subsampling_lengths(2), N/subsampling_lengths(3)];

% find subsampled DFTs
f_1 = fft(f(1:subsampling_rates(1):end));
f_2 = fft(f(1:subsampling_rates(2):end));
f_3 = fft(f(1:subsampling_rates(3):end));

% debug
[B_1, ~] = sort(abs(f_1), 'descend');
[B_2, ~] = sort(abs(f_2), 'descend');
[B_3, ~] = sort(abs(f_3), 'descend');

%threshold for noise robustness
thresh1 = 0.5 * max(abs(f_1));
thresh2 = 0.5 * max(abs(f_2));
thresh3 = 0.5 * max(abs(f_3));

%find peaks
peaks1 = find(abs(f_1) > thresh1);
peaks2 = find(abs(f_2) > thresh2);
peaks3 = find(abs(f_3) > thresh3);

x = crt([peaks1(1)-1, peaks2(1)-1, peaks3(1)-1],subsampling_lengths); 
%convert from 1- to 0- index for CRT above

 disp(x) %revert to 0-index