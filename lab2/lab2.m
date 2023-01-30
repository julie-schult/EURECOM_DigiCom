%% 2: OFDM demodulation
%extract the SSS symbol based on the timing information obtained from the PSS. 
%Remove the prefix and perform the FFTs of size 2048 to transfer them to the frequency domain.
%{
r3 = exp(-j*2*pi*Mmax*(0:length(rxs3)-1)*delta_f/f_s)'.*rxs3;

rpss = r3(N_f+144:2048);
rsss = r3(N_f+2192+2192+(1:2048));

% Removing prefix
rm_pref = r3(145:length(r3));

%FFT
Rpss = fft(rm_pref);
Rsss = fftshift(fft(rsss));
%}

rpss = rxs0(N_f+144+(0:length(pss2_t)-1));
rsss = rxs0(N_f+2*2192+144+(0:2048-1));
rsss1 = rxs1(N_f+2*2192+144+(0:2048-1));
rsss2 = rxs2(N_f+2*2192+144+(0:2048-1));
rsss3 = rxs3(N_f+2*2192+144+(0:2048-1));
Rpss = (fft(rpss,2048));
Rsss = fftshift(fft(rsss,2048));
Rsss1 = fftshift(fft(rsss1,2048));
Rsss2 = fftshift(fft(rsss2,2048));
Rsss3 = fftshift(fft(rsss3,2048));
%% 2a: Plot the complex modulus of the 127 resource elements corresponding to the SSS.
a=find(Rsss>10000);
Rsss_new = Rsss(a(2):a(length(a)));
figure(1)
stem(abs(Rsss_new)), title('Rsss')
Rsss0 = Rsss_new(56:182);
figure(2)
stem(abs(Rsss0)), title('Rsss')
%Characteristic of QAM OFDM: ?

%% 2a Rsss1
figure()
stem(abs(Rsss1))
%% 2a Rsss2
a = find(Rsss2>10000);
Rsss2_new = Rsss2(a(2):a(length(a)));
figure()
stem(abs(Rsss2_new)), title('Rsss2')
Rsss22 = Rsss2_new(61:187);
%% 2a Rsss3
figure()
stem(abs(Rsss3))
a = find(Rsss3>10000);
Rsss3_new = Rsss3(a(2):a(length(a)));
figure()
stem(abs(Rsss3_new)), title('Rsss3')
Rsss33 = Rsss3_new(73:199);
%% 2b
figure(3)
plot(real(Rsss0),imag(Rsss0),'bo') 
%The 'circle' of points is the sss part of the OFDM symbol

%% 2c
%figure(4)
%stem(abs(fftshift(Rpss)))

b=find(Rpss>4000); %???
Rpss_new = Rpss(961:1087);%???

H = conj(pss_2).*Rpss;
H = [H(1985:2048); H(1:63)];
figure(4)
plot(abs(H)), title('freq domain channel est.')
figure(5)
h = ifft(H);
%v = [-1024:1023];
plot(abs(h)), title('time domain channel est.') %Strong peak in 0, this means this is a good channel estimator

%% 2d
%{
Mlr = 2*real(H.*Rpss_new.*fft(H.*Rpss_new));

figure(6)
plot(Mlr)
%}
%% 2e


%% 2f
iq_sss0 = conj(H).*Rsss0;
iq_sss2 = conj(H).*Rsss22;
iq_sss3 = conj(H).*Rsss33;

figure(7)
plot(real(iq_sss0),imag(iq_sss0),'o'), title('I/Q signal Rsss0')
figure(8)
plot(real(iq_sss2),imag(iq_sss2),'o'), title('I/Q signal Rsss2')
figure(9)
plot(real(iq_sss3),imag(iq_sss3),'o'), title('I/Q signal Rsss3')


%% 3
res = d(1+2:3:size(d,1),:)*iq_sss2;
[val,id] = max(normalize(abs(res)))
figure(10)
stem(abs(res)), title('Correlation BPSK ans SSS0')

%% 4
cellID = 2 + 3*id
