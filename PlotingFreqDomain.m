%% duo Ca plot thru ffy on freq. domain
Fs = 5000;            % Sampling frequency = 5kHz
TR = 1/Fs;
nTR = round(length(x1)/Fs); %numeber of TR will be a approximate value into integer
load rsfmri_13.mat
x1 = channels{1,8}.data;
x2 = channels{1,9}.data;
t = 0:1/Fs:(length(x1)-1)*TR;
L1 = numel(channels{1,3}.data);% Length of signal
L2 = numel(channels{1,4}.data);
% Y1 = fft(channels{1,3}.data);
% Y2 = fft(channels{1,4}.data);
% NFFT1 = 2^nextpow2(L1);% Next power of 2 from length of y           
% NFFT2 = 2^nextpow2(L2);
% Y1_fft = fft(x1,NFFT1)/L1;
% Y2_fft = fft(x2,NFFT2)/L2;
% x1 = 2*abs(Y1_fft(1:NFFT1/2+1));
% x2 = 2*abs(Y2_fft(1:NFFT2/2+1));
% f1 = Fs/2*linspace(0,1,NFFT1/2+1);
% f2 = Fs/2*linspace(0,1,NFFT2/2+1);
x1_mean = x1 - mean(x1);
x2_mean = x2 - mean(x2);
% hold on;
% plot(f1,2*abs(Y1_fft(1:NFFT1/2+1)));
% plot(f2,2*abs(Y2_fft(1:NFFT2/2+1)));
% grid on;
% hold off;
% xlabel('frequency(Hz)');
% ylabel('Amplitude');
plot(t,x1_mean,t,x2_mean);% plot averaged data of 2 filtered channel
xlim([0 60]);
xlabel('Time(s)');
ylabel('Amplitude');
title('Time domain signal expression');

norm_x1 = zscore(x1_mean);
norm_x2 = zscore(x2_mean);

plot(t,norm_x1,t,norm_x2);
% plot(t,x1_mean,t,x2_mean);
xlim([0 60]);
xlabel('Time(s)');
ylabel('Amplitude');
title('Time domain signal expression');
legend('ca_01','ca_02');

[Pxx1,w1] = pwelch(norm_x1,[],[],[],Fs);
[Pxx2,w2] = pwelch(norm_x2,[],[],[],Fs);
plot(w1,10*log10(Pxx1),w2,10*log10(Pxx2));
xlim([0 10]);
xlabel('freq.(Hz)')
ylabel('Power(dB)')
% plot(w1,pow2db(Pxx1),w2,pow2db(Pxx2));% same plot as 10*log10()
legend('Ca_01','Ca_02');

%% Filtered plot showing identical features
windowSize = 70; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;
y_1 = filter(b,a,10*log10(Pxx1));
y_2 = filter(b,a,10*log10(Pxx2));
plot(w1,y_1,w2,y_2);
xlim([0 10]);
xlabel('freq.(Hz)')
ylabel('Power(dB)')
title('Smoothened plot')
legend('ca01','ca02')

% BPF with .01-.1 Hz:
Fc_low = 1;
Fc_high = 6;
% N = 46;
N = 60;
beta = 0.005;
win = kaiser(N+1, beta);
flag = 'scale';  % Sampling Flag

% Calculate the BPF-ed pow2freq  using the FIR1 function.
b  = fir1(N, [Fc_low, Fc_high]/(Fs/2), 'bandpass', win, flag);

fvtool(b,1,'Fs',Fs);
% Uncomment for freq. domain:(2 lines)
% bpf1 = filter(b,1,y_1);
% bpf2 = filter(b,1,y_2);
bpf1 = filter(b,1,x1_mean);
bpf2 = filter(b,1,x2_mean);
% Uncomment for freq. domain:(1 line)
% plot(w1,bpf1,w2,bpf2)
plot(t,bpf1,t,bpf2)
% Uncomment for freq. domain:(5 lines)
% xlim([0 10]);
% xlabel('freq.(Hz)')
% ylabel('Power(dB)')
% legend('ca01','ca02');
% title('1-6Hz BPF 2channels pow2freq');
xlim([0 60]);
xlabel('Time(sec)')
ylabel('Amplitude')
legend('ca01','ca02');
title('1-6Hz BPF 2channels timecourse');


[Pxx1,w1] = pwelch(bpf1,[],[],[],Fs);
[Pxx2,w2] = pwelch(bpf2,[],[],[],Fs);
plot(w1,10*log10(Pxx1),w2,10*log10(Pxx2));
xlim([0 10]);
xlabel('freq.(Hz)')
ylabel('Power(dB)')
% hold on
% plot(w1,pow2db(Pxx1),w2,pow2db(Pxx2));% same plot as 10*log10()
legend('Ca_01','Ca_02');
title('pow2freq_1_6BPF');

%% not necessary below
freq_1 = meanfreq(Pxx1,w1);% mean freq. estimate of channels
freq_2 = meanfreq(Pxx2,w2);
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% f = Fs*(0:(L/2))/L;
% plot(Y1) ;
% title('Single-Sided Amplitude Spectrum of X(t)')
% xlabel('f (Hz)')
% ylabel('|P1(f)|')

                    
%% PSD way

pspectrum(x1_mean,Fs,'spectrogram','FrequencyResolution',0.5)% here 0.5Hz is resolution of frequency axis
% ylim([0 10]);
window = hann(2^13);
ovrl = 1/2*length(window);
nfft = 2^18;
[s,f,t] = spectrogram(abs(x1_mean),window,ovrl,nfft,Fs,'yaxis');
imagesc(t(t>0),f(f<10),abs(s(f<10)));
ax = gcs;
ax.YDir = 'reverse';
xlabel('Time(s)');
ylabel('Freq.(kHz)');
colorbar;
colormap(jet);

%% Attemp w/ 3rd-party filter
ca_matrix = x1_mean;
windowSize = 2;
% tDuration = TR*nTR;
tDuration = 473;


[tf freqs times] = timefreq(ca_matrix,Fs,'wletmethod','dftfilt3','winsize',windowSize*Fs,'ntimesout',round(tDuration/windowSize*2),'freqs',[0,20]);
%              tf      = complex time frequency array for all trials (freqs, times, trials)
%              freqs   = vector of computed frequencies (Hz)
%              times   = vector of computed time points (ms)
    
  figure;      % subplot (311);
        tf = abs(tf); % get the power
        times = times/1000; %mse --> sec
        imagesc(times,freqs,tf);
%         set(gca,'ydir','normal','xlim',[times(1) times(end)],'ylim',[0.1 20]);
%         xlabel(sprintf('Time in sec')); 
        ylabel('Frequency in Hz'); grid on;
        set(gca,'XTick',[0:60:60*tDuration]);
%         set(gca,'XTickLabel',[0:60:60*tDuration]);
        title(strcat('calcium spectrum'));
        varname = 'raw_ca';
        %title(strcat(num2str(minNum(1,iMin)),'min after surgery'));  %%%%%
        %for automation use the iMin
        h=colorbar;
        ylabel(h,'Power');
        caxis([0, 30]); 
        grid off;
        title(strcat(varname,'after PBF()'));
        set(gcf, 'PaperUnits', 'inches');
        y_width=7 ;x_width=20;
        set(gcf, 'PaperPosition', [0 0 x_width y_width]); % 
        saveas(gcf,strcat(varname,'-0p1to20Hz',name),'tiff'); 
        close all;