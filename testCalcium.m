load('C:\Users\wangqi\Documents\Lab\Data\ca_pupil_08212018\rsfmri_13.mat')
addpath("C:\Users\wangqi\Downloads\Matlab_scripts\");
Fs = 5000;
TR = 1/Fs;
Spatial_Resolution = 0.05;
%% channel01 and 02 calcium data
raw_ca_01 = channels{1,3}.data;
raw_ca_02 = channels{1,4}.data;
x1 = raw_ca_01(1,:);
x2 = raw_ca_02(1,:);
a = 1;
b = [1/4 1/4 1/4 1/4];
y1 = filter(b,a,x1);%moving average filter
y2 = filter(b,a,x2);
t = 1:length(x1);

t_1 = 0:1/Fs:(length(x1)-1)/Fs;
plot(t_1,x1);
xlabel('Time (seconds)');
ylabel('Amplitude');
title('Raw data timecourse');
%% plot power 2 freq.
m_1 = length(x1);       % original sample length
n_p1 = pow2(nextpow2(m_1));  % transform length, exponential 2 is to facilitate calculation
fft_x1_p = fft(x1,n_p1);
f_p1 = (0:n_p1-1)*(Fs/n_p1)/10;% DFT of signal
power = abs(fft_x1_p).^2/n_p1;   
pwr2db = 10*log10(power);
figure,
plot(f_p1(1:floor(n_p1/2)),power(1:floor(n_p1/2)));% how to adjust intensity range of power for display??
xlim([0 5]);
ylim([0 5]);
xlabel('Frequency');
ylabel('Power');

subplot(2,1,1);
plot(t,x1,'--',t,y1,'-');
xlim([0 length(x1)*TR]);
ylim([-1.5,-1.0]);
title("Raw data Calcium Channel 01 VS Averaged Calcium data");
legend('Original Data','Locally Averaged Data');
subplot(2,1,2);
plot(t,x2,'--',t,y2,'-');
xlim([0 length(x1)*TR]);
title("Raw data Calcium Channel 02 VS Averaged Calcium data");
legend('Original Data','Locally Averaged Data');
% kspace_raw(:,:) = fft2c(raw_ca(:,:));
% imshow(abs(kspace_raw(:,:)));
plot(t,y1,'--',t,y2,'-'); % plot averaged data
xlim([0 4751022*TR]);
title('Illustration of filtered calcium data from both channels')
legend('Channel01 filetred','Channel02 filtered');
% averaged channels calcium
x1_ave = x1-mean(x1);
x2_ave = x2-mean(x2);
figure,subplot(2,1,1);
plot(t,x1_ave);
xlim([0 150]);
xlabel('time(t)');
ylabel('averaged amplitude')
subplot(2,1,2);
plot(t,x2_ave);
xlim([0 150]);
xlabel('time(t)');
ylabel('averaged amplitude');
%% pwelch method to draw freq. domain figure
x = x1;
[pxx,f] = pwelch(x,[],[],[],Fs, 'onesided'); 
h= figure();
% plot(f,10*log10(pxx),'LineWidth',plot_linewidth);
plot(f,pow2db(pxx))
%% Straight forward way to freq. domain
n_x1 = length(x1_ave);% number of samples
fft_x1 = fft(x1_ave);
f = (0:n_x1-1)*(Fs/n_x1);     % frequency range
power = (abs(fft_x1).^2)/n_x1;    % power of the DFT

plot(f,power)
xlabel('Frequency')
ylabel('Power')
xlim([0 5]);
ylim([-10 10]);
%% PSD 
window = hann(2^13);
ovrl = 1/2*length(window);
nfft = 2^18;
[s,f,t] = spectrogram(x1_ave,window,ovrl,nfft,Fs);

psd_x1 = imagesc(t,f(f<10),abs(s(f<10)));
colormap(jet);
arr_psd_x1 = imread(psd_x1);
Col_psd_x1(:,1) = mean(arr_psd_x1(:,1));

