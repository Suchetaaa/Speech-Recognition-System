[s, fs] = audioread('s8.wav');
hanningWindow = hanning(256);
nfft = 256;
noverlap = 156;
figure;
[S,F,T] = spectrogram(s, hanningWindow, noverlap, nfft, fs);
imagesc(T,F,log(abs(S)))
set(gca,'YDir','Normal')
xlabel('Time (secs)')
ylabel('Freq (Hz)')
title('Short-time Fourier Transform spectrum')