[s,fs] = audioread('s8.wav');
sound(s, fs);
[f0, loc] = pitch(s,fs);    %f0 has the pitch coefficients s
plot(loc, f0);