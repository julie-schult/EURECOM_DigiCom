% this is a 5G PRACH format 0 transmission
L=839;

% this is number of samples per bin, floor(L/Ncs) gives the number of cyclic shifts, see below
Ncs=26;
% this is the FFT size for the generation/reception of PRACH
N=49152;
% this is the length of the cyclic prefix for PRACH
Ncp=6336;


% 6-bit data messages for 3 transmitters / UEs
preamble_index1=63;
preamble_index2=31;
preamble_index3=11;

% up to 6 Zadoff-Chu root sequences for this format
utab=[129 710 140 699 120 719];
% number of cyclic shifts
nshifts = floor(L/Ncs)
% number of Zadoff-Chu sequences required
nseq = ceil(64/nshifts)

% index of the preamble sequence to use
uind1=floor(preamble_index1/nshifts)
uind2=floor(preamble_index2/nshifts)
uind3=floor(preamble_index3/nshifts)
% index of cyclic shift to use
nuind1=rem(preamble_index1,nshifts)
nuind2=rem(preamble_index2,nshifts)
nuind3=rem(preamble_index3,nshifts)
if (uind1>=length(utab) || uind2>=length(utab) || uind3>=length(utab)) 
  fprintf("ERROR tab length %d : %d %d %d",length(utab),uind1,uind1,uind3) 
end

% These are the Zadoff-Chu Sequence generators (time-domain) 
% for the 3 transmitters
xu1 = exp(-j*pi*utab(1+uind1)*(0:838).*(1:839)/839);
xu2 = exp(-j*pi*utab(1+uind2)*(0:838).*(1:839)/839);
xu3 = exp(-j*pi*utab(1+uind3)*(0:838).*(1:839)/839);

% implement cyclic-shifts
% Note we do this in the time-domain and then do an 839-point fft here in MATLAB
% This is not usually done in practice because of the complexity of the FFT (i.e. a large prime number)
% There is a way to compute the Fourier transform directly and then perform the cyclic shift by a multiplication of a phasor in the frequency-domain.

xuv1 = xu1; 
xuv2 = xu2;
xuv3 = xu3;
for (n=0:838)
  xuv1(n+1) = xu1(1+rem(n+(Ncs*nuind1),839));
  yuv1 = fft(xuv1);
  xuv2(n+1) = xu2(1+rem(n+(Ncs*nuind2),839));
  yuv2 = fft(xuv2);
  xuv3(n+1) = xu3(1+rem(n+(Ncs*nuind3),839));
  yuv3 = fft(xuv3);
end

% put the PRACH in the lowest frequency (positive) subcarriers starting at carrier 7
Xuv1 = zeros(1,49152);
Xuv1(7+(1:839)) = yuv1;
Xuv2 = zeros(1,49152);
Xuv2(7+(1:839)) = yuv2;
Xuv3 = zeros(1,49152);
Xuv3(7+(1:839)) = yuv3;

% bring to time-domain
xuv1_49152 = ifft(Xuv1);
xuv2_49152 = ifft(Xuv2);
xuv3_49152 = ifft(Xuv3);

% add cyclic prefix
xuv1_49152 = [xuv1_49152((49152-6335):end) xuv1_49152];
xuv2_49152 = [xuv2_49152((49152-6335):end) xuv2_49152];
xuv3_49152 = [xuv3_49152((49152-6335):end) xuv3_49152];

% normalizes the transmit signal to unit-energy
xuv1_49152 = xuv1_49152/sqrt(sum(abs(xuv1_49152).^2)/length(xuv1_49152));
en1=mean(abs(xuv1_49152).^2)
xuv2_49152 = xuv2_49152/sqrt(sum(abs(xuv2_49152).^2)/length(xuv2_49152));
en2=mean(abs(xuv2_49152).^2)
xuv3_49152 = xuv3_49152/sqrt(sum(abs(xuv3_49152).^2)/length(xuv3_49152));
en3=mean(abs(xuv3_49152).^2)

% Plot the time-domain and frequency-domain waveform (xuv1)
figure(1)
plot(10*log10(abs(xuv1_49152))), title('Time-domain waveform')

figure(2)
plot(10*log10(abs(fftshift(fft(xuv1_49152))))), title('Frequency-domain waveform')

% Question: What can you say regarding the frequency span (approximately how many PRBs does this waveform occupy
%The frequemcy span is calculated below in Hz
span = length(find(10*log10(abs(fftshift(fft(xuv1_49152))))>25));
%The number of PRBs are given by
size_PRB = 12*30; % since mu is 1 we have 12 subcarriers * 30
PRBs = round(span/size_PRB); % 3

% simulate time-delay
delay1=300;
delay2=140;
delay3=40;
delaymax = 1+max([delay1 delay2 delay3]);
xuv1_49152 = [zeros(1,delay1) xuv1_49152 zeros(1,delaymax-delay1)];
xuv2_49152 = [zeros(1,delay2) xuv2_49152 zeros(1,delaymax-delay2)];
xuv3_49152 = [zeros(1,delay3) xuv3_49152 zeros(1,delaymax-delay3)];

SNR=0;
snr=10.^(.1*SNR);
noise1 = sqrt(.5/snr)*(randn(1,length(xuv1_49152))+sqrt(-1)*randn(1,length(xuv1_49152)));
noise2 = sqrt(.5/snr)*(randn(1,length(xuv1_49152))+sqrt(-1)*randn(1,length(xuv1_49152)));
rxsig1_justnoise = xuv1_49152 + noise1;
rxsig2_justnoise = xuv1_49152 + xuv2_49152 + xuv3_49152 + noise2;

% do TDL-C channel generation
fs=61.44e6;
SCS=30e3;
DS=300e-9;

H=get_tdl(fs,SCS,[0:105],DS,'tdlc');
H2 = zeros(1,2048);
halflength =53*12;
H2((2048-(halflength-1)):2048) = H(1:halflength);
H2(1:halflength) = H(halflength+(1:halflength));
h1 = ifft(H2)*sqrt(2048);

H=get_tdl(fs,SCS,[0:105],DS,'tdlc');
H2 = zeros(1,2048);
halflength =53*12;
H2((2048-(halflength-1)):2048) = H(1:halflength);
H2(1:halflength) = H(halflength+(1:halflength));
h2 = ifft(H2)*sqrt(2048);

H=get_tdl(fs,SCS,[0:105],DS,'tdlc');
H2 = zeros(1,2048);
halflength =53*12;
H2((2048-(halflength-1)):2048) = H(1:halflength);
H2(1:halflength) = H(halflength+(1:halflength));
h3 = ifft(H2)*sqrt(2048);

rxsig3_noiseandchannel = conv(h1,xuv1_49152);
rxsig3_noiseandchannel = rxsig3_noiseandchannel + sqrt(.5/snr)*(randn(1,length(rxsig3_noiseandchannel))+sqrt(-1)*randn(1,length(rxsig3_noiseandchannel)));

rxsig4_noiseandchannel = conv(h1,xuv1_49152) + conv(h2,xuv2_49152) + conv(h3,xuv3_49152);
rxsig4_noiseandchannel = rxsig4_noiseandchannel + sqrt(.5/snr)*(randn(1,length(rxsig4_noiseandchannel))+sqrt(-1)*randn(1,length(rxsig4_noiseandchannel)));



% What to do now
% a) implement the receiver using a frequency-domain correlation using the Zadoff-Chu sequences generation method as above

%https://www.rfwireless-world.com/5G/5G-NR-Zadoff-chu-sequence.html
%https://se.mathworks.com/help/lte/ug/synchronization-signals-pss-and-sss.html
X_uv = [ exp(-j*pi*utab(1)*(0:L-1).*(1:L)/L); exp(-j*pi*utab(2)*(0:L-1).*(1:L)/L); ]; %Channel, utab(1) and utab(2) because nseq=2

%% Rxsig1
%close all
rxsig1_remove_prefix = rxsig1_justnoise(Ncp + 1 : N); %Removing prefix (Ncp=6336)
rxsig1_fft = fft(rxsig1_remove_prefix,N); %N-point fft
rxsig1_freq_cut = rxsig1_fft((7+(1:L))); % k bar is 7

X1_uv_freq = fft(X_uv(1,:), L); %First column of X_uv (with utab(1))
X2_uv_freq = fft(X_uv(2,:), L); %Second column of X_uv (with utab(2))

rxsig1_freq_corr1 = rxsig1_freq_cut.*conj(X1_uv_freq); %Correlated in frequency domain with the first column of the channel
rxsig1_freq_corr2 = rxsig1_freq_cut.*conj(X2_uv_freq); %Correlated in frequency domain with the second column of the channel

rxsig1_time1 = ifft(rxsig1_freq_corr1,L); % inverse fourier transform
rxsig1_time2 = ifft(rxsig1_freq_corr2,L); % inverse fourier transform

figure(3)
plot(abs(rxsig1_time1));
%title("Correlation of rxsig1");
hold on;
plot(abs(rxsig1_time2));
title("Correlation of rxsig1 with the receiver");
xlabel("Lag") 
ylabel("Correlation") 
legend("Correlation 1", "Correlation 2");
hold off;

%Explanation: We can see from the figure that the correlation 2 has a strong peak a
%little after 0 (lag 40). Here, the lag is 0, which tells us that the rxsig1 is received by the 2nd
%Zadoff-chu sequence. Nothing for correlation 1.

%% Rxsig2
%close all
rxsig2_remove_prefix = rxsig2_justnoise(Ncp + 1 : N); %Removing prefix (Ncp=6336)
rxsig2_fft = fft(rxsig2_remove_prefix,N); %N-point fft
rxsig2_freq_cut = rxsig2_fft((7+(1:L))); % k bar is 7

rxsig2_freq_corr1 = rxsig2_freq_cut.*conj(X1_uv_freq); %Correlated in frequency domain with the first column of the channel
rxsig2_freq_corr2 = rxsig2_freq_cut.*conj(X2_uv_freq); %Correlated in frequency domain with the second column of the channel

rxsig2_time1 = ifft(rxsig2_freq_corr1,L); % inverse fourier transform
rxsig2_time2 = ifft(rxsig2_freq_corr2,L); % inverse fourier transform

figure(4)
plot(abs(rxsig2_time1));
hold on;
plot(abs(rxsig2_time2));
title("Correlation of rxsig2 with the receiver");
xlabel("Lag") 
ylabel("Correlation") 
legend("Correlation 1", "Correlation 2");

%Explanation: We can see from the figure that the correlation 1 has a strong peak a
%little after 0 (lag 40)  and again at approx. lag 555 which tells us that the rxsig2 is received by the 1st
%Zadoff-chu sequence. Correlation 2 also has a peak a little after 0 which
%means rxsig2 is also received by the 2nd column. (?????)

%% Rxsig3
%close all
rxsig3_remove_prefix = rxsig3_noiseandchannel(Ncp + 1 : N); %Removing prefix (Ncp=6336)
rxsig3_fft = fft(rxsig3_remove_prefix,N); %N-point fft
rxsig3_freq_cut = rxsig3_fft((7+(1:L))); % k bar is 7

rxsig3_freq_corr1 = rxsig3_freq_cut.*conj(X1_uv_freq); %Correlated in frequency domain with the first column of the channel
rxsig3_freq_corr2 = rxsig3_freq_cut.*conj(X2_uv_freq); %Correlated in frequency domain with the second column of the channel

rxsig3_time1 = ifft(rxsig3_freq_corr1,L); % inverse fourier transform
rxsig3_time2 = ifft(rxsig3_freq_corr2,L); % inverse fourier transform

figure(5)
plot(abs(rxsig3_time1));
hold on;
plot(abs(rxsig3_time2));
title("Correlation of rxsig3 with the receiver");
xlabel("Lag") 
ylabel("Correlation") 
legend("Correlation 1", "Correlation 2");

%Similar to rxsig1, but stronger peak for correlation 2 and a little more noise for
%correlation 1, but no peak.

%% Rxsig4
%close all
rxsig4_remove_prefix = rxsig4_noiseandchannel(Ncp + 1 : N); %Removing prefix (Ncp=6336)
rxsig4_fft = fft(rxsig4_remove_prefix,N); %N-point fft
rxsig4_freq_cut = rxsig4_fft((7+(1:L))); % k bar is 7

rxsig4_freq_corr1 = rxsig4_freq_cut.*conj(X1_uv_freq); %Correlated in frequency domain with the first column of the channel
rxsig4_freq_corr2 = rxsig4_freq_cut.*conj(X2_uv_freq); %Correlated in frequency domain with the second column of the channel

rxsig4_time1 = ifft(rxsig4_freq_corr1,L); % inverse fourier transform
rxsig4_time2 = ifft(rxsig4_freq_corr2,L); % inverse fourier transform

figure(6)
plot(abs(rxsig4_time1));
hold on;
plot(abs(rxsig4_time2));
title("Correlation of rxsig4 with the receiver");
xlabel("Lag") 
ylabel("Correlation") 
legend("Correlation 1", "Correlation 2");

%Similar to rxsig2, but now, correlation 1 has a way stronger peak at lag 555,
%and the peak at lag 40 is the same.
%For correlation 2 there is only one peak at lag 40, and this correlation
%is stronger than for rxsig2.

% b) show how the data detection and time-delay estimation
%% Time-delay estimation of correlations (prob. not needed for every correlation?)

%Time-delay of rxsig1 correlation 1
[~, lag] = max(abs(rxsig1_time1)) %index (lag) and value of the maximum peak in correlation
bin = floor(lag/nshifts); %binary index of peak
sig_start = (bin-1)*Ncs - 1;
delay = lag - sig_start;
time_delays(1) = delay;


%Time-delay of rxsig2 correlation 1
[~, lag] = max(abs(rxsig2_time1))
bin = floor(lag/nshifts);
sig_start = (bin-1)*Ncs - 1;
delay = lag - sig_start;
time_delays(2) = delay;

%Time-delay of rxsig3 correlation 1
[~, lag] = max(abs(rxsig3_time1))
bin = floor(lag/nshifts);
sig_start = (bin-1)*Ncs - 1;
delay = lag - sig_start;
time_delays(3) = delay

%Time-delay of rxsig4 correlation 1
[~, lag] = max(abs(rxsig4_time1));
bin = floor(lag/nshifts);
sig_start = (bin-1)*Ncs - 1;
delay = lag - sig_start;
time_delays(4) = delay

%Time-delay of rxsig1 correlation 2
[~, lag] = max(abs(rxsig1_time2));
bin = floor(lag/nshifts);
sig_start = (bin-1)*Ncs - 1;
delay = lag - sig_start;
time_delays(5) = delay;

%Time-delay of rxsig2 correlation 2
[~, lag] = max(abs(rxsig2_time2));
bin = floor(lag/nshifts);
sig_start = (bin-1)*Ncs - 1;
delay = lag - sig_start;
time_delays(6) = delay;

%Time-delay of rxsig3 correlation 2
[~, lag] = max(abs(rxsig3_time2));
bin = floor(lag/nshifts);
sig_start = (bin-1)*Ncs - 1;
delay = lag - sig_start;
time_delays(7) = delay;

%Time-delay of rxsig4 correlation 2
[~, lag] = max(abs(rxsig4_time2));
bin = floor(lag/nshifts);
sig_start = (bin-1)*Ncs - 1;
delay = lag - sig_start;
time_delays(8) = delay

