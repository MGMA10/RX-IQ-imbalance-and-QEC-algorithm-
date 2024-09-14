%% Variables Initialization
clear;
carrier_frequancy = 1e3; % Carrier frequency in Hz
input_frequancy = 10; % Input signal frequency in Hz
sambels = 0.1 * [1+i,1-i,-1-i,-1+i]; % Modulation symbols
gain_imbalance = 1.122; % Gain imbalance for demodulation
phase_imbalance = 0; % Phase imbalance in radians

% Select a random symbol from the constellation
samb = sambels(randi(4));

% Define time step for carrier and input signals
time = 0:1/(carrier_frequancy)/100:1; 

% Generate cosine and sine waves for carrier and input signals
cos_range_fc = cos(2*pi*carrier_frequancy*time);
sine_range_fm = cos(2*pi*input_frequancy*time);

figure;
plot(time, sin(2*pi*input_frequancy*time) +cos(2*pi*input_frequancy*time));
title('Signal in time Domain');
xlabel('time (s)');
ylabel('Magnitude');



%% Modulation

    % Modulate signal (real part on cosine, imaginary part on sine)
    cos_branch_mod = real(samb) * cos(2*pi*input_frequancy*time) .* cos(2*pi*carrier_frequancy*time);
    sin_branch_mod = imag(samb) * -1 * sin(2*pi*input_frequancy*time) .* sin(2*pi*carrier_frequancy*time);
  
mod_signal = cos_branch_mod + i * sin_branch_mod;
N = length(mod_signal);




% Perform FFT and plot the frequency domain representation
fft_mod_signal = fftshift(abs(fft(mod_signal)));

figure;
plot(-50000:50000, fft_mod_signal);
title('Modulated Signal in Frequency Domain');
xlabel('Frequency (Hz)');
ylabel('Magnitude');


%% Channel (Optional, We can add noise or distortion here)

%% Demodulation
cos_branch_demod = mod_signal .* cos(2*pi*carrier_frequancy*time) * 2;
sin_branch_demod = gain_imbalance * 2 * mod_signal .* sin(2*pi*carrier_frequancy*time + phase_imbalance);


out_dem = cos_branch_demod + j * sin_branch_demod ;

gain_bar = rms(sin_branch_demod)/rms(cos_branch_demod);

phase_bar = sum( cos_branch_demod .* sin_branch_demod)/sqrt(sum(sin_branch_demod .^2)* sum(cos_branch_demod .^2));
phase_bar1=asin(phase_bar);

% Perform FFT and plot the frequency domain representation
fft_mod_signal = fftshift(abs(fft(cos_branch_demod)));

figure;
plot(-50000:50000, fft_mod_signal);
title('Modulated Signal in Frequency Domain');
xlabel('Frequency (Hz)');
ylabel('Magnitude');


out = lowpass(out_dem,10,carrier_frequancy);
h=fir1(1024,10/carrier_frequancy/2);
h_len = length(h);
con=conv(out_dem,h);
outfft=fftshift(abs(fft(con)));

% Plot demodulated signal in the frequency domain
figure;
plot (-50000:50000, fftshift(abs(fft(out))));
title('demodulated signal in the frequency domain');
xlabel('time (s)');
ylabel('Magnitude');


% Plot demodulated signal in time domain
figure;
plot(abs(con(floor(h_len/2):N-floor(h_len/2))));
title('Input and Demodulated Signal in Time Domain');
xlabel('Time (s)');
ylabel('Amplitude');
% hold on;
% plot(sin(2*pi*input_frequancy*time) + cos(2*pi*input_frequancy*time));
% plot(sin(2*pi*input_frequancy*time));

v=fftshift(abs(fft(out)));

IQ_Imbalance = -20*log10(max( v(51000:53000 ))/(max( v(47000:49000))))

cos_branch_demod_new = cos_branch_demod;

sin_branch_demod_new = tan(phase_bar1) .* cos_branch_demod + sin_branch_demod ./ (gain_bar * cos(phase_bar1));

out_corr = cos_branch_demod_new + j * sin_branch_demod_new;
% Plot demodulated signal in the frequency domain after correction
figure;
plot (-50000:50000, fftshift(abs(fft(out_corr))));
title('demodulated signal in the frequency domain after correction');
xlabel('time (s)');
ylabel('Magnitude');

w=fftshift(abs(fft(out_corr)));

IQ_Imbalance_corr = -20*log10(max( w(51000:53000 ))/(max( w(47000:49000))))
% nfft = 2.^ nextpow2 (N);
% fy = fft (out_dem , nfft);
% fy = fy (1:nfft/2);
% xfft = carrier_frequancy/2 .* (0:nfft/2-1)/nfft;
% 
% figure
% subplot(2,2,1);
% plot(xfft,(fy /max(fy)));
% subplot(2,2,2);
% plot(time, out_dem);
% h=fir1(128,15/carrier_frequancy/2);
% con=conv(out_dem,h);
% subplot(2,2,3);
% 
% hold on;
% plot(con);