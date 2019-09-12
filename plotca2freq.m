%% I/O
loc = strcat('C:\Users\wangqi\Documents\Lab\Data\ca_Hang');
numb = 17;
% for num=numb(1:length(numb));
num = numb;
cd(loc);
FileName_data  = strcat('Scan',num2str(num),'_TK','.mat');
% uncomment for Patricia's data
% FileName_data  = strcat('rsfmri_',num2str(num),'.mat');
load(FileName_data);
%% duo Ca plot thru ffy on freq. domain
linewidth = 1;
Fs = channels{1,7}.samples_per_second;            % Sampling frequency 
TR = 1/Fs;
% load rsfmri_13.mat
x1 = channels{1,7}.data;
nTR = round(length(x1)/Fs); %numeber of TR will be a approximate value into integer
% load rsfmri_13.mat
% x2 = channels{1,9}.data;
t = 0:1/Fs:(length(x1)-1)*TR;
L1 = numel(channels{1,3}.data);% Length of signal
% L2 = numel(channels{1,4}.data);
x1_mean = x1 - mean(x1);
var_1 = zscore(lowpass(x1_mean,20,Fs));
var_2 = normalize(lowpass(x1_mean,20,Fs),'range',[0 1]);
% var_1 = x1_mean;
% stimuli = channels{1,4}.data;
% sti_fs = channels{1,4}.samples_per_second;
% stu_t = 0:1/sti_fs:(length(stimuli)-1)/sti_fs;
% x2_mean = x2 - mean(x2);
% plot(t,x1_mean,t,x2_mean);% plot averaged data of 2 filtered channel
maxsec_ca  = floor(length(var_1)/Fs);
% maxsec_sti = floor(length(stimuli)/sti_fs);
% plot(t,x1_mean)
hold on
subplot(2,1,1),
plot(var_1(1:maxsec_ca*Fs),'linewidth',linewidth);
xlim([1 maxsec_ca*Fs]);
%ylim([1500 8000]);
set(gca,'XTick',[0:20*Fs:maxsec_ca*Fs]);
set(gca,'XTickLabel',[0:20:maxsec_ca]);
% xlim([0 60]);
xlabel('Time(s)');
ylabel('Amplitude(<20Hz)');
title('Time domain normalized signal');
grid on;
subplot(2,1,2),
plot(var_2(1:maxsec_ca*Fs),'linewidth',linewidth);
% plot(stimuli(1:maxsec_sti*sti_fs/2),'linewidth',linewidth);
xlim([1 maxsec_ca*Fs]);
%ylim([1500 8000]);
set(gca,'XTick',[0:20*Fs:maxsec_ca*Fs]);
set(gca,'XTickLabel',[0:20:maxsec_ca]);
% xlim([0 60]);
xlabel('Time(s)');
ylabel('Amplitude(<20Hz)');
title('Time domain averaged signal');
legend('ca','stimuli')
grid on;
%% freq. domain w/ smoothen plot
[Pxx1,w1] = pwelch(var_1,[],[],[],Fs);
windowSize = 5; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;
y = filter(b,a,Pxx1);
plot(w1,10*log10(y));% beside 6Hz is a cardiac spike
xlim([0 10]);
xlabel('Frequency(Hz)')
ylabel('Power(dB)');
title('Power to Frequency below 10Hz')
% plot(w1,pow2db(Pxx1),w2,pow2db(Pxx2));% same plot as 10*log10()
% legend('Ca_01','Ca_02');

%% Filtered plot showing identical features
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
bpf1 = filter(b,1,y);
% bpf2 = filter(b,1,y_2);
%uncomment for time domain:(1 line)
% bpf1 = filter(b,1,x1_mean);
% bpf2 = filter(b,1,x2_mean);

% Uncomment for freq. domain:(1 line)
% plot(w1,bpf1,w2,bpf2)
bpf1 = mean(bpf1);
plot(w1,bpf1);
% plot(t,bpf1,t,bpf2)
%uncomment for time domain:(1 line)
% plot(t,bpf1)
% Uncomment for freq. domain:(5 lines)
xlim([0 10]);
xlabel('freq.(Hz)')
ylabel('Power(dB)')
% legend('ca01','ca02');
title('1-6Hz BPF  pow2freq');
% uncomment for time domain:(5 lines)
% xlim([0 60]);
% xlabel('Time(sec)')
% ylabel('Amplitude')
% legend('ca01','ca02');
% title('1-6Hz BPF 2channels timecourse');


%% not necessary below
freq_1 = meanfreq(Pxx1,w1);% mean freq. estimate of channels
% freq_2 = meanfreq(Pxx2,w2);

% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% f = Fs*(0:(L/2))/L;
% plot(Y1) ;
% title('Single-Sided Amplitude Spectrum of X(t)')
% xlabel('f (Hz)')
% ylabel('|P1(f)|')

                    
%% PSD way

[P,F,T] = pspectrum(var_2,Fs,'spectrogram','FrequencyResolution',0.1,'FrequencyLimits',[0.1 10]);% here 0.5Hz is resolution of frequency axis
% p = abs(P);
% imagesc(T,F,P);
% pow_min = min(P(:));
% % [rrr,ccc] = min(P(:));
% % [rrr,ccc] = ind2sub(size(P),ccc);
% pow_max = max(P(:));
% set('ylim',[0.1 10]);
% caxis([0,0.1]);
% colorbar

% pspectrum(var_2,Fs,'spectrogram','FrequencyResolution',0.1,'FrequencyLimits',[0.1 10]);
% how to tune color span to fit current range w/ rgb colors ???




window = hann(2^12);
ovrl = 1/2*length(window);
nfft = 2^18;
[s,f,t] = spectrogram(abs(x1_mean),window,ovrl,nfft,Fs,'yaxis');
imagesc(t(t>0),f(f<10),abs(s(f<10)));
ylim([0 10]);
xlabel('Time(s)');
ylabel('Freq.(kHz)');
colorbar;
colormap(jet);

%% Attemp w/ 3rd-party filter
% ca_matrix = var_1;
ca_matrix = x1_mean;
windowSize = 2;
% tDuration = TR*nTR;
tDuration = maxsec_ca;
addpath('/Users/wangqi/Downloads/Matlab_scripts/calcium/');

[tf freqs times] = timefreq(ca_matrix,Fs,'wletmethod','dftfilt3','winsize',windowSize*Fs,'ntimesout',round(tDuration/windowSize*2),'freqs',[0,20]);
%              tf      = complex time frequency array for all trials (freqs, times, trials)
%              freqs   = vector of computed frequencies (Hz)
%              times   = vector of computed time points (ms)
    
  figure;      % subplot (311);
        tf = abs(tf); % get the power
        times = times/1000; %mse --> sec
        imagesc(times,freqs,tf);
        ylabel('Frequency(Hz)'); grid on;
        xlabel('Time(sec)');
        set(gca,'ydir','normal','xlim',[times(1) times(end)],'ylim',[0.1 20]);% reversed y_axis
        %         set(gca,'XTick',[0:60:60*tDuration]);
        title('Calcium spectrum(Task-Related)');
        varname = 'raw_ca';
        h=colorbar;
        ylabel(h,'Power(dB)');
        caxis([0, 10]); 
        grid off;
        set(gcf, 'PaperUnits', 'inches');
        y_width=7 ;x_width=20;
        set(gcf, 'PaperPosition', [0 0 x_width y_width]); %  
        
%         norm_tf_temp = abs(tf(:,:))./max(abs(tf(:,:)));
%         Threshold_tf = abs(tf(:,:))./mean(norm_tf_temp(:,:),2);
%         norm_tf = abs(Threshold_tf(:,:))./max(max(abs(Threshold_tf(:,:))));

%% seperated bands by frequency span
        [Ix0p5 C0p5]   = find(freqs' == 0.5);
        [Ix1p0 C1p0]   = find(freqs' == 1);
        [Ix2p0 C2p0]   = find(freqs' == 2);
        [Ix3p0 C3p0]   = find(freqs' == 3);
        [Ix4p0 C4p0]   = find(freqs' == 4);
        [Ix5p0 C5p0]   = find(freqs' == 5);
        [Ix6p0 C6p0]   = find(freqs' == 6);
        [Ix7p0 C7p0]   = find(freqs' == 7);
        [Ix8p0 C8p0]   = find(freqs' == 8);
        [Ix9p0 C9p0]   = find(freqs' == 9);
        [Ix10p0 C10p0] = find(freqs' == 10); 
       
        Avg_0p5to1= mean(tf(Ix0p5:Ix1p0,:)); %band1
        Avg_1to2  = mean(tf(Ix1p0:Ix2p0,:)); %band2
        Avg_2to3  = mean(tf(Ix2p0:Ix3p0,:)); %band3
        Avg_3to4  = mean(tf(Ix3p0:Ix4p0,:)); %band4
        Avg_4to5  = mean(tf(Ix4p0:Ix5p0,:)); %band5
        Avg_5to6  = mean(tf(Ix5p0:Ix6p0,:)); %band6
        Avg_6to7  = mean(tf(Ix6p0:Ix7p0,:)); %band7
        Avg_7to8  = mean(tf(Ix7p0:Ix8p0,:)); %band8
        Avg_8to9  = mean(tf(Ix8p0:Ix9p0,:)); %band9
        Avg_9to10 = mean(tf(Ix9p0:Ix10p0,:)); %band10

        Avg_freqs_decom = [Avg_0p5to1', Avg_1to2', Avg_2to3', Avg_3to4', Avg_4to5', Avg_5to6', Avg_6to7', Avg_7to8', Avg_8to9', Avg_9to10']; 
        MostAvg_freqs_decom = [Avg_0p5to1', Avg_1to2', Avg_2to3', Avg_5to6', Avg_6to7', Avg_7to8', Avg_8to9', Avg_9to10']; 
        [m,n]=size(MostAvg_freqs_decom);
        for i=1:n
            [envMostAvg_freqs_decom(:,i),del] = envelope(MostAvg_freqs_decom(:,i));
        end
        ln=1;
        
        figure;
                subplot(611)
        plot(Avg_freqs_decom(:,1),'linewidth',ln);ylabel('0.5-1');xlim([1 m]);set(gca,'xticklabel',[]);
        title('frequency in seperated span');
                subplot(612)
        plot(Avg_freqs_decom(:,2),'linewidth',ln);ylabel('1-2');xlim([1 m]);set(gca,'xticklabel',[]);
                subplot(613)
        plot(Avg_freqs_decom(:,3),'linewidth',ln);ylabel('2-3');xlim([1 m]);set(gca,'xticklabel',[]);
                subplot(614)
        plot(Avg_freqs_decom(:,7),'linewidth',ln);ylabel('6-7');xlim([1 m]);set(gca,'xticklabel',[]);
                subplot(615)
        plot(Avg_freqs_decom(:,9),'linewidth',ln);ylabel('8-9');xlim([1 m]);set(gca,'xticklabel',[]);
                subplot(616)
        plot(Avg_freqs_decom(:,10),'linewidth',ln);ylabel('9-10');xlim([1 m]);
% end       
%% all span plot normalized[Revise later]
% [rw,col] = find(tf(8:16,1:30) >= 2);
start_row = find(abs(freqs-3)<0.0001); % the signal of the row to plot
start_col = 20;
% windowSize = 5; 
% b = (1/windowSize)*ones(1,windowSize);
% a = 1;
% tf_filtered = filter(b,a,tf(start_row,start_col:t_num));
% tf = tf - mean(tf);
var_3 = tf ;% this variable is matrix to plot
for t_num = start_col:666

if mod(t_num,20) == 0
%     minu_t = t_num - start_col
    hold on;
    grid on;
    plot(var_3(start_row,start_col:start_col+20))
    xlim([1 20]);
    xlabel('time(sec)');
    ylabel('Power(dB)');
    start_col = start_col+20;
else
disp('Not a period')    
end
end
% story of 3Hz vs 6Hz
plot(tf(start_row,:));
xlabel('Times(sec)');
ylabel('Amplitude');
xlim([0 maxsec_ca]);
title('Averaged calcium signal in 3Hz');
hold on;
plot(tf(start_row_6,:));
xlabel('Times(sec)');
ylabel('Amplitude');
xlim([0 maxsec_ca]);
title('Averaged calcium signal in 3Hz');
legend('3Hz response','6Hz Harmonic');
[pk_1,lc_1] = findpeaks(tf(start_row,:),'MinPeakHeight',5);
plot(lc_1,pk_1,'r','LineWidth',2);
xlabel('Times(sec)');
ylabel('Amplitude');
title('Peak value in line');
% try 6Hz for comparison
start_row_6 = find(abs(freqs-6)<0.0001);
% timecourse of 6Hz
[pk_2,lc_2] = findpeaks(tf(start_row_6,:),'MinPeakHeight',2);
plot(lc_2,pk_2,'b','LineWidth',2);
xlabel('Times(sec)');
ylabel('Amplitude');
title('Peak value in line');


% averaged PSD on time
 ave_psd = mean(tf(:,times(1):times(end)),2);
 plot(ave_psd);
 xlabel('Frequency(Hz)');
 ylabel('PSD(dB)');
 xticks(0:4:80);
 xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20'})
 yticklabels({'0','1','2','3','4','5','6','7','8','9','10'});
 title('Averaged PSD');
 grid on;
%% convolve
%  %%
%  numb = [20 31];
%  for i = 1:length(numb)
% Ca_Band(:,i)=dataCa(:,band,i); % store a new variable with only the desired band, for each run. Units are in averaged power (power acquired with timefreq).
% end
%  h = colorbar;
% oldsr=1; % from seconds..
% newsr=1/1.5; % ..to 1/TR
% [P,Q] = rat(newsr/oldsr);
% rs_h=resample(h,P,Q); % figure,plot(rs_h);
% 
% for i=1:length(numb)
% newvector(:,i)=conv(rs_h,Ca_Band(:,i));
% end
% Ca_Band_conv=newvector(1:length(Ca_Band),:);
%     figure,
%     subplot(211)
%     plot(Ca_Band(:,1), 'b');
%     legend({'band'});
%     subplot(212)
%     plot(Ca_Band_conv(:,1), 'r');
%     legend({'convolved band'});