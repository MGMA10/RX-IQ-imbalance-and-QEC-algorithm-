clear;
f = 10 ;
gain_imbalance = 1.1; % Gain imbalance for demodulation
phase_imbalance = 0.06; % Phase imbalance in radians
sam_f = 1e3;

% Define time step for carrier and input signals
time = 0:1/(sam_f):1; 

I = cos(2*pi*f*time);
Q = gain_imbalance * sin (2*pi*f*time - phase_imbalance);

out = I + j * Q ;

w=fftshift(abs(fft(out)));
IQ_Imbalance = -20*log10(max( w(1:500))/(max( w(500:1000))))

%%
% Perform FFT and plot the frequency domain representation
fft_mod_signal = fftshift(abs(fft(out)));

%% correction 
gain_bar = rms(Q)/rms(I);

phase_bar = -sum( I .* Q)/sqrt(sum(Q .^2)* sum(I .^2));
% phase_bar1=asin(phase_bar);

I_new = I;

Q_new = tan(phase_bar) * I + Q / (gain_bar * cos(phase_bar));

out_corr = I_new + j * Q_new;

v=fftshift(abs(fft(out_corr)));
IQ_Imbalance_corr = -20*log10(max( v(1:500))/(max( v(500:1000))))

%%
% Perform FFT and plot the frequency domain representation
fft_mod_signal_new = fftshift(abs(fft(out_corr)));
figure;
plot(-500:500, abs(fft_mod_signal_new));
title('Modulated Signal in Frequency Domain After Correction');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([-40 40])

figure
plot(time, imag(out_corr));
title('Signal in time Domain After Correction');
hold on
plot(time, imag(out));
xlabel('Time');
ylabel('Magnitude');
legend('After Correction' ,'Before Correction' );

figure;
plot(-500:500, fft_mod_signal);
title('Modulated Signal in Frequency Domain Before Correction');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([-40 40])