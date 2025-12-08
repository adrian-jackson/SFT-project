%% SFT Project:
% Implement the SFT algorithm in MATLAB. Make comparisons with the built-in 
% FFT algorithm of MATLAB in terms of running time as the length of data 
% increases, for different sparsity cases

% enable the following to show plots, etc., but will slow the program down
debug_flag = true;

% set vars
filter_taps = 16; %must be less than signal length
signal_length = 70;
amount_of_peaks = 2;
center_frequency_peak_index = 20;
%% create signal
% need to be able to make signals with N = length,k = number of peaks, 
% p = spacing of peaks, epsilon[j] = noise function

K = amount_of_peaks; 
A = 1;
% use the folloing 3 lines for gaussian-like magnitudes
%Amax = K; % peak magnitude for central sample
%x = 1:K;
%A = K/Amax + (Amax - K/Amax) * (1 - abs(x - (K+1)/2)/((K+1)/2-1));
k0 = center_frequency_peak_index; % center frequency index
N = signal_length; % signal length
delta_w = 0.6; % peak spacing, (0 to 1]. approach 0 for closer spacing
n  = 0:N-1; % time indices

k_list = (1:N/K:N) + ((k0-512) + N/(2*K));
half_k_range = round(((k_list(end) - k_list(1)) * delta_w) / 2);   % half of total range
k_list = linspace(k0 - half_k_range, k0 + half_k_range, K);
k_list = round(k_list);
k = zeros(1,N);
k(k_list) = A;

%plot freq peaks
if debug_flag
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
end
% Create time signal
S = exp(2*pi*1i * (k_list(:) * n) / N);
f = sum(S, 1);    % 1-by-N signal

%% Implement SFT I. Filter to Isolate Frequencies

filter_shifts = 3;
f_filtered_bank = zeros(K,N,filter_shifts);

for l = 1:1:filter_shifts
    fprintf('Shifting filters by exp(1j * 2*pi*1/(%d*2K) \n', l-1)
    M = filter_taps; % number of FIR taps
    b_fir = fir1(M-1, 1/(K*2)); % simple FIR LP
    h = b_fir(:).'; 
    hn = 1:1:M;
    shift = N/K; %need to shift K filter copies by N/K each
    %f_filtered = zeros(K,N); unnecessary for toeplitz
    h_bank = zeros(K, M);
    for shifter = 0:1:K-1
        freq_shift = pi*shifter/K;
        h_bank(shifter+1,:) = h .* exp(1j * freq_shift * hn) .* exp(1j * 2*pi*1/(K*4) * hn) .* exp(1j * pi*((l-1)/(4*K)) * hn);
    end
    
    if debug_flag
        figure; hold on;
        cmap = lines(K);  % generate K distinct colors
        for j = 1:K
            [H, w] = freqz(h_bank(j,:), 1, 2048);  % frequency response of j-th filter
            plot(w/pi, abs(H), 'LineWidth', 1.5, 'Color', cmap(j,:));
        end
        xlabel('Normalized Frequency (\times \pi rad/sample)');
        ylabel('Magnitude');
        title('Shifted FIR Filterbank Responses');
        grid on;
        legend(arrayfun(@(k) sprintf('Shift %d', k-1), 1:K, 'UniformOutput', false));
    end
    
    % multiply K filters of length M=128 by 1xN signal f using Toeplitz matrix
    % Build Toeplitz matrix for convolution
    F = zeros(N,M);
    for m = 1:M
        F(m:end, m) = f(1:end-m+1);
    end
    
    % Multiply h_bank by X' (rows = filters)
    f_filtered = h_bank * F.';  % K x N
    f_filtered_bank(:,:,l) = f_filtered;
    
    % debug check: did binning isolate one frequency per subband signal?
    if debug_flag
        for j=1:1:K
            F_temp = fft(f_filtered(j,:));
            [~,peak_temp] = max(abs(F_temp));
            fprintf('Peaks found by filtering in subband %d: ', j);
            fprintf('%.2f ', peak_temp-1);
            fprintf('\n');       
        end
    end
end


%% Implement SFT II. Identify Frequency Peaks

% Aliasing and CRT
% calculate padded signal of length with specific coprime factors.
% subsampling rates MUST be > k but << N

subsampling_rates = [2,5,7];
remainders = zeros(1,length(subsampling_rates));
for l=1:1:filter_shifts
    for s=1:1:K  %NEED TO EXLUDE FILTERS THAT FIND >1 PEAK
        fc = f_filtered_bank(s,:,l);

        for i = 1:length(subsampling_rates)
            nk = subsampling_rates(i);
            
            % Subsample the signal every nk samples
            F_sub = fft(fc(1:(N/nk):end));
                      
            % Find peak in subsampled FFT
            [~, idx] = max(abs(F_sub));
            
            % Convert MATLAB 1-based index to 0-based remainder
            remainders(i) = idx-1;
            
            fprintf('Subsample %d: peak at index %d (mod %d)\n', i, remainders(i), nk);
        end
        k_recovered = crt(remainders, subsampling_rates);
        fprintf('Recovered frequency index using CRT: %d\n', k_recovered); %matlab 1-indexing 
    end
end