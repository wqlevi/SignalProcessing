load rsfmri_13.mat

fs = channels{1,3}.samples_per_second; % sampling rate
tr = 1/fs;
mkdir ./Processing_Test/MATLAB/test_result_pupile
cd ./Processing_Test/MATLAB/test_result_pupile  
figure();
raw_ca = plot(channels{1,3}.data);% raw data
title('Raw Calcium Data Time course');
xlabel('Time course numbering');
ylabel('Intensity');
saveas(raw_ca,'Raw Calcium Data Time Domain','png');

ca = channels{1,3}.data - mean(channels{1,3}.data); % set baseline to be 0
figure();
bsl_ca = plot(ca);
title('0-Baseline Calcium Data Time course');
xlabel('Time(s)');
ylabel('Intensity');
saveas(bsl_ca,'0-Shifted Calcium Data Time domain','png');

[Pxx1,w1] = pwelch(channels{1,3}.data,[],[],[],fs); % do pwelch to plot a PSD
figure();
subplot(2,1,1);
freq_raw_ca = plot(w1,pow2db(Pxx1)); % transfer Pxx1 to power unit (dB)
title('Raw Calcium Data Frequency domain');
xlabel('Frequency(Hz)');
ylabel('Power(dB)');
grid on;
subplot(2,1,2);
freq_raw_ca = plot(w1,pow2db(Pxx1)); % transfer Pxx1 to power unit (dB)
xlim([0 10]);
title('Raw Calcium Data Frequency domain');
xlabel('Frequency(Hz)');
ylabel('Power(dB)');
grid on;
saveas(freq_raw_ca,'Raw Calcium Data Frequency domain','png');

[Pxx1,w1] = pwelch(ca,[],[],[],fs); % do pwelch on extracted data
figure();
subplot(2,1,1);
grid on;
freq_shi_ca = plot(w1,pow2db(Pxx1));
title('0-Shifted Calcium Data Frequency domain');
xlabel('Frequency(Hz)');
ylabel('Power(dB)')
subplot(2,1,2);
freq_shi_ca = plot(w1,pow2db(Pxx1));
title('0-Shifted Calcium Data Frequency domain');
xlim([0 10]);
xlabel('Frequency(Hz)');
ylabel('Power(dB)')
saveas(freq_shi_ca,'0-Shifted Calcium Data Frequency domain','png');
% This function below showed raw freq. plot :
% figure, plot(w1,Pxx1)

% figure(), plot(w1,10*log10(Pxx1))% same thing as pow2db(Pxx1)
% title('Calcium Data Frequency domain');
% xlabel('Frequency(Hz)');
% ylabel('Power(dB)')
%% setting for PSD windows
win = hann(2^13); % hanning window size, the larger the smaller resolution
ovrl = 1/2 * length(win);
nfft = 2^18; % resolution of pwelch result horizontally
%% Plot PSD
[S,F,T] = spectrogram(ca,win,ovrl,nfft,fs);% what're S,F,T for?
figure
psd_ca_power = imagesc(T, F(F<150), abs(S((F<150), :)));
title('Calcium Data Frequency domain');
xlabel('Time Course(s)');
yticks(0:1:50);
ylabel('Frequncy(Hz)')
ylim([0 4]);
c = colorbar;
c.Label.String = 'Power Intensity(dB)';
set(gcf,'color','w');
colormap(jet);
saveas(psd_ca_power,'Calcium Data Frequency domain','png');

figure;
psd_ca_log = imagesc(T, F(F<150), 10*log10(abs(S((F<150), :))));
title('Calcium Data Frequency domain');
xlabel('Time Course(s)');
ylabel('Frequncy(Hz)')
ylim([0 4]);
c = colorbar;
c.Label.String = 'Power Intensity(dB)';
saveas(psd_ca_log,'Calcium Data Frequency domain','png');
