Nfft=2048;
pss_0;
pss_1;
pss_2;

pss0_t = ifft(pss_0.'); 
pss0_t = pss0_t/norm(pss0_t);

pss1_t = ifft(pss_1.');
pss1_t = pss1_t / norm(pss1_t);

pss2_t = ifft(pss_2.');
pss2_t = pss2_t / norm(pss2_t);

pss0_t = [pss0_t(((2048-143):2048)) pss0_t];
pss1_t = [pss1_t(((2048-143):2048)) pss1_t];
pss2_t = [pss2_t(((2048-143):2048)) pss2_t];

%% 1.1.1 Question 1: plot of real and imaginary part, and magnitude
%Choose to work on pss0
figure(1)
subplot(3,1,1)
plot(real(pss0_t)), title('Real part of pss0')
subplot(3,1,2)
plot(imag(pss0_t)), title('Imaginary part of pss0')
subplot(3,1,3)
plot(20*log10(abs(pss0_t))), title('Magnitude of pss0')

%% 1.1.1 Question 2: Power spectrum and bandwidth
%Continue working on pss0
Fs=61.44*10.^6;
len=2192;
Sx=fftshift(20*log10(abs(fft(pss0_t))));
l=length(find(Sx>-10));
BW=(l*Fs)/len;

figure(2)
plot(fftshift(20*log10(abs(fft(pss0_t))))), title('Power spectrum on a dB scale')
%{
We estimated the bandwidth to be around 4 MHz
%}

%% 1.1.1 Question 3: Auto- and crosscorrelation of the three signals
[Auto1, lag1]=xcorr(pss0_t);
[Auto2, lag2]=xcorr(pss1_t);
[Auto3, lag3]=xcorr(pss2_t);

[cross01, lag01]=xcorr(pss0_t,pss1_t);
[cross02, lag02]=xcorr(pss0_t,pss2_t);
[cross12, lag12]=xcorr(pss1_t,pss2_t);

figure(3)
subplot(321)
plot(lag1,real(Auto1)), title('Real part of autocorrelation of pss0')
subplot(322)
plot(lag1,imag(Auto1)), title('Imaginary part of autocorrelation of pss0')
subplot(323)
plot(lag2,real(Auto2)), title('Real part of autocorrelation of pss1')
subplot(324)
plot(lag2,imag(Auto2)), title('Imaginary part of autocorrelation of pss1')
subplot(325)
plot(lag3,real(Auto3)), title('Real part of autocorrelation of pss2')
subplot(326)
plot(lag3,imag(Auto3)), title('Imaginary part of autocorrelation of pss2')

figure(4)
subplot(321)
plot(lag01,real(cross01)), title('Real part of crosscorrelation of pss0 and pss1')
subplot(322)
plot(lag01,imag(cross01)), title('Imaginary part of crosscorrelation of pss0 and pss1')
subplot(323)
plot(lag02,real(cross02)), title('Real part of crosscorrelation of pss0 and pss2')
subplot(324)
plot(lag02,imag(cross02)), title('Imaginary part of crosscorrelation of pss0 and pss2')
subplot(325)
plot(lag12,real(cross12)), title('Real part of crosscorrelation of pss1 and pss2')
subplot(326)
plot(lag12,imag(cross12)), title('Imaginary part of crosscorrelation of pss1 and pss2')

%{
Since the amplitude around 0 is almost 0 for all crosscorrelations, we can
say that the signals are almost orthogonal (If the amplitude is exactly 0
at 0 the signals are orthogonal).

If we assume that they are orthogonal, the ratio signal energy to
interference in dB will be the signal energy??

%}

%% 1.2.1 Question 2, Time and frequency plots on the dB scale
load('rxsignal_justnoise.mat');
load('rxsignal_withchannel.mat');
load('rxsignal_withchannelandfreqoff.mat');

figure(5)
subplot(311)
plot(20*log10(abs(rxs0))), title('time representation of rxs0')
subplot(312)
plot(20*log10(abs(rxs1))), title('time representation of rxs1')
subplot(313)
plot(20*log10(abs(rxs2))), title('time representation of rxs2')

figure(6)
subplot(311)
plot(fftshift(20*log10(abs(fft(rxs0))))), title('frequency representation of rxs0')
subplot(312)
plot(fftshift(20*log10(abs(fft(rxs1))))), title('frequency representation of rxs1')
subplot(313)
plot(fftshift(20*log10(abs(fft(rxs2))))), title('frequency representation of rxs2')

%% 1.2.1 Question 3
% Choosing to look at rxs2
Fs=61.44*10.^6;
len=614400;
Sx0=fftshift(20*log10(abs(fft(rxs2))));
l0=length(find(Sx0>100));
BW0=(l0*Fs)/len;

%{
The bandwidth is estimated to be 7.2 MHz
Outside the band of interest there are noise and interference components.

1.2.1 Question 4
The signal changes shape when its affected by the channel.
%}

%% 1.3.1 

%{

magnitude square of the duration of the signal N_pss samples (2192)
x*(n)y(n) - projection of received signal on the transmit signal
y(n) on x(n), this is a convolution of x* and y(-n)|n=N_f

non-coherent receiver, because we have an unknown channel (like unknown phase in the continuous case)
y(n + N_f) delaying received signal by N_f (time) to get the correlation between
transmit signal and received signal, e^(delta_f) is the frequency offset

The expression is basically a correlation, in time and frequency, and we
want to maximize the correlation to receive the best possible signal

plot(abs(conv(x_i,conj(fliplr(y)))))
-> should see a strong peak in N_f

Question 1

Question 2

%}

%% 1.3.1 Question 3


figure(7)
subplot(311)
plot(abs(conv(fliplr(conj(pss0_t)),rxs3)))
subplot(312)
plot(abs(conv(fliplr(conj(pss1_t)),rxs3)))
subplot(313)
plot(abs(conv(fliplr(conj(pss2_t)),rxs3)))

%We see a strong peak for pss2, which means that it is most likely the transmittet signal
%Finding N_f
[max_value,max_pos] = max(abs(conv(fliplr(conj(pss2_t)),rxs3)));
length_pss2 = length(pss2_t);
N_f = max_pos - (length_pss2 - 1);


%% 1.3.1 Question 4

r_pss = rxs3(N_f+(0:length(pss2_t)-1));
delta_f = 100;
f_s = 61.44e6;

statistic = zeros(151,2192);
%n = 0:length(r_pss)-1:

for m = -75:75
    statistic(m+76,:) = conj(pss2_t).*exp(-2*pi*i*(0:length(r_pss)-1)*m*(delta_f/f_s));
end    

Y = abs(statistic*r_pss).^2;

[max_value,max_pos] = max(Y)
Mmax = max_pos-75; 
off = Mmax*delta_f

figure(8)
plot(Y);


