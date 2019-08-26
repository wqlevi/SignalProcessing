load('C:\Users\wangqi\Documents\Lab\Data\ca_pupil_08212018\rsfmri_16.mat')
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
%% PSD 
window = hann(2^13);
ovrl = 1/2*length(window);
nfft = 2^18;
[s,f,t] = spectrogram(bpf_x1,window,ovrl,nfft,Fs);

imagesc(t,f(f<10),abs(s(f<10)));
colormap(jet);

