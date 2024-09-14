clear;
f = 10 ;
gain_imbalance = 1.122; % Gain imbalance for demodulation
phase_imbalance = 0; % Phase imbalance in radians
sam_f = 1e3;

% Define time step for carrier and input signals
time = 0:1/(sam_f):1; 
i=0;
for gain_imbalance = logspace(-3,3,10000)
    i=i+1;
I = cos(2*pi*f*time);
Q = gain_imbalance * sin (2*pi*f*time - phase_imbalance);

out = I + j * Q ;

w=fftshift(abs(fft(out)));
IQ_Imbalance(i) = 20*log10(max( w(1:500))/(max( w(500:1000))));
end

plot(20*log10(logspace(-3,3,10000)),IQ_Imbalance)
title('IQ Imbalance gain')
xlabel('gain_imbalance , dB');
ylabel('Image Rejection ,dB');
grid on
xlim([-3 3])
ylim([-50 -10])


i=0;
gain_imbalance = 1;
for phase_imbalance = -0.2:0.001:0.2
    i=i+1;
I = cos(2*pi*f*time);
Q = gain_imbalance * sin (2*pi*f*time - phase_imbalance);

out = I + j * Q ;

w=fftshift(abs(fft(out)));
IQ_Imbalancep(i) = 20*log10(max( w(1:500))/(max( w(500:1000))));
end
figure;
plot((-0.2:0.001:0.2)*180/pi,IQ_Imbalancep)
title('IQ Imbalance phase')
xlabel('Phase_imbalance , dB');
ylabel('Image Rejection ,dB');
grid on
xlim([-10 10])
 ylim([-50 -20])
 
