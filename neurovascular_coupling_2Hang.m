%%%%%%%%%%%%%%%%%%%%%% USAGE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script aims to operate the following functions:
% - frequency dependent cross-correlation between calcium and fmri;(<0.05Hz & 2-3Hz Ca)
% - Statistical plot of BOLD along cortical depth;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fmri_path = 'D:\05232019\05232019_fmri_analyzed';
% fmri_path = 'D:\05302019_Hang\05302019_analyzed';
fmri_path = 'D:\09122019\fmri_analyzed09122019';
% fmri_path = 'D:\10062018_evoke\1006fmri\10062018';
% fmri_path = 'D:\03242019\fmri_analyzed03242019';
% fmri_path='D:\HANG\10062018rsesting-state\fmri_data';
% fmri_path = 'D:\10092018\fmri_analyzed10092018';% 10092018
% ca_path = 'D:\05232019\05232019_ca';
% ca_path = 'D:\05302019_Hang\05302019_ca';
ca_path = 'D:\09122019\09122019_ca';
% ca_path = 'D:\10062018_evoke\ca_1006\0_Task_related';
% ca_path = 'D:\03242019\Calcium_data03242019';
% ca_path = 'D:\HANG\10062018rsesting-state\calcium_data';
% ca_path = 'D:\10092018\ca_10092018';% 10092018
mkdir('C:\Users\wangqi\Desktop\Neurovascular coupling\0912\figs');
savepath = 'C:\Users\wangqi\Desktop\Neurovascular coupling\0912\figs';
addpath('C:\Users\wangqi\Documents\Lab\Demo\Calcium');

TR = 0.1;% 0.1s
Spatial_res = 0.05; %0.05 per pixel/voxel
prestim = 10;%10 sec
fmri_duration = 640;
bsl_fmri = 2;%sec
% scans = [35 37 39 44 49 51 57 59 64 66 69];%03242019
% scans = [22,26];
% scans = [15 18 20 23 25 28 31 34 35 37 39 41 43 45 48 49 52 54 56 58 59];%05232019TK
% scans = [36,51,55,57,59,61,63,65,69,71];%10062018
% scans = [16 17 20 23 25 27 29 31 33 35 36 ];%05302019tk_P
% scans = [13,15,17,19,21,23,25,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64,66,68,70,72,74,76,78];%09122019_RS
scans = [12,14,16,18,20,22,24,27,29,31,33,35,37,39,41,43,45,47,49,51,53,55,57,59,61,63,65,67,69,71,73,75,77,79];%09122019_TK
% scans = 12;
% scans = [39 41 43 45 47 49 51 53 55 57 59 61 63 65 67 69 71 73 75 77 79 81];%05302019tk_N
% scans = [16 17 20 23 25 27 29 31 33 35 36 39 41 43 45 47 49 51 53 55 57 59 61 63 65 67 69 71 73 75 77 79 81];%05302019tk
% scans = [15 19 22 24 26 28 30 32 34 37 38 40 42 44 46 48 50 52 54 56 58 60 62 64 66 68 70 72 74 76 78 80];%05302019
% scans = [17 19 21 22 24 26 27 32 33 36 38 40 42 44 46 48 50 51 53 55 57];%05232019
% scans = [24,29,31,35,37,39,41,44,46,48,50,52,54,56,58,60,62,64,66,68,70,72,74];%10062018_TK
fmri_dummy = ones(fmri_duration/TR,1);
chan_ca = 7;
nslice = 3;
nvoxel = 2/Spatial_res;
cortical_depth = 2/Spatial_res;
offset = 1;
evoke = 1;
nfreq = 5;
voxel_idx = [1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5];
%% load signal

% load calium
for ir = 1:length(scans)
    data_ca_orig = load([ca_path,'\scan_',num2str(scans(ir)),'.mat']);
    [data_match,fmri_dummy,beg,fin] = match_acq_fmri(data_ca_orig,fmri_dummy,TR,prestim);
    ca_match{ir} = -(data_match.channels{chan_ca}.data)';
%     baseline_ca{ir} = mean(-data_ca_orig.channels{7}.data(:,1:prestim*fs_ca),2); 
end
% set time scales
fs_ca =  data_match.channels{chan_ca}.samples_per_second;
t_ca = [0:1/fs_ca:(length(ca_match{ir})-1)/fs_ca];
t_fmri = [0:TR:(length(fmri_dummy)-1)*TR];ls
% load fmri
for ir = 1:length(scans)
    for is = 1:nslice
        tmp = load([fmri_path, '\', num2str(scans(ir)), '\Results_Slice', num2str(is), ...
                '\line_scanning_data_s', num2str(is), '.mat']);
         fmri{ir}{is} = abs(tmp.total_cortical_depth_map);   
    end
end
%% reshaping pre-process
cd (savepath);

for ir = 1:length(scans)
    for is = 1:nslice
        mean_slice = mean(fmri{ir}{is}, 2);
        mean_slice = mean_slice - min(mean_slice);% mean of one scan
        
        half_intens = max(mean_slice)/2;
        i_half(is,ir) = find(mean_slice>half_intens, 1);% the first voxel with value larger than trail average
        
        % shift the cortex lower by several voxels
        i_half(is,ir) = i_half(is,ir)+offset;
    end
end
% trim 2mm of cortex
for ir = 1:length(scans)
    for is = 1:nslice
        fmri_trim{ir}{is} = fmri{ir}{is}(i_half(is,ir)+1:i_half(is,ir)+cortical_depth,:);
    end
end

% tSNR & cortical map
for ir = 1:length(scans)
    for is = 1:nslice
        h = figure;
        imagesc(t_fmri,1:size(fmri{ir}{is},1),fmri{ir}{is});
        hold on;
        colormap jet;
        xlabel('time(s)');
        ylabel('cortical depth');
        yline(i_half(is,ir),'b--','LineWidth',2);
        yline(i_half(is,ir)+cortical_depth,'b--','LineWidth',2);
        title(['Cortical map trail ',num2str(ir),' slice ',num2str(is)]);
    end
end
for ir = 1:length(scans)
    for is = 1:nslice
        h = figure;
        plot(mean(fmri{ir}{is},2));
        xline(i_half(is,ir),'b--','LineWidth',2);
        xline(i_half(is,ir)+cortical_depth,'b--','LineWidth',2);
        xlim([0,128]);
        title(['SNR trail ',num2str(ir),' slice ',num2str(is)]);
    end 
end

%detrending ca & fmri
for ir = 1:length(scans)
    for is = 1:nslice
        fmri_tc{ir}{is} = detrend(mean(fmri_trim{ir}{is},1));
    end
    ca_match{ir} = detrend(ca_match{ir});
end

%% Raw signals 1-D plot
% fmri time course each slice
for ir = 1:length(scans)
    for is = 1:nslice
        h = figure;
%         fmri_tc_nor{ir}{is} = (fmri_tc{ir}{is} - mean(fmri_tc{ir}{is}(:,1:bsl_fmri/TR),2))/abs(mean(fmri_tc{ir}{is}(:,1:bsl_fmri/TR),2));
        fmri_tc_nor{ir}{is} = zscore(fmri_tc{ir}{is});
        plot(t_fmri,fmri_tc_nor{ir}{is},'r-');
        xlim([0,t_fmri(end)]);
        ylim([-5,5]);
        xlabel('time(s)');
        ylabel('fMRI a.u.');
        set(gcf,'Position',[500 500 980 300 ]);
        title(['fmri trail ',num2str(scans(ir)),' slice ',num2str(is)]);
        box off;
        saveas(h,['fmri raw trail ',num2str(scans(ir)),' slice ',num2str(is),'.jpg']);
    end
end
% demean ca itself
for ir = 1:length(scans)
        ca_demean{ir} =  ca_match{ir} - mean(ca_match{ir});
        ca_per{ir} = (ca_match{ir} - mean(ca_match{ir}))/abs(mean(ca_match{ir}));
%         ca_demean{ir} = (ca_match{ir} - baseline_ca{ir})/baseline_ca{ir};
        h = figure;
%         plot(t_ca,zscore(ca_demean{ir}));
        plot(t_ca,ca_per{ir});
        xlim([0,t_ca(end)]);
        xlabel('time(s)');
        ylabel('Ca^2^+ \Delta F/F');
        set(gcf,'Position',[500 500 980 300 ]);
        title(['Ca^2^+ timecourse trail ',num2str(scans(ir))]);
        saveas(h,['Ca timecourse trail ',num2str(scans(ir)),'.jpg']);
end
close all
%% spectral analysis(spectrogram & timefreq)
% ca spectrogram
disp('psd ca...')
win = 2*fs_ca;%2sec
NFFT = 2^18;
overlap = 0;
len_ca = length(ca_demean{ir});

for ir = 1:length(scans)
    for c = 1:4
        [tf1(:,:,c),freq1,time1(:,c)] = timefreq(ca_demean{1}((len_ca/4)*(c-1)+1:c*(len_ca/4),:),fs_ca,'wletmethod','dftfilt3','winsize',win,'ntimesout',6400/4,'padratio',256,'freqs',[0,10]);
    end
    save(['tf',num2str(ir),'.mat'], 'tf1');
    clear tf1
end

for ir = 1:length(scans)
    h = figure;
tf2{ir} = reshape(tf1{ir},size(freq1,2),6400);
freq2 = reshape(freq1,[],1);
time2 = reshape(time1,[],1)/1000;
imagesc(time2,freq2,abs(tf2{ir}));
colormap jet;
caxis([0,6]);
set(gca,'YDir','normal');
xlabel('time(s)');
ylabel('Frequency(Hz)');
ylim([freq2(1),freq2(end)]);
set(gcf,'Position',[500 500 980 300 ]);
saveas(h,['ca 0.01-10Hz PSD trail ',num2str(scans(ir)),'.jpg']);

ca_pp_bsl{ir} = abs(tf2(freq2<0.05,:));% baseline assignment(<0.05Hz)
    ca_pp_bls_ave{ir} = mean(ca_pp_bsl{ir},1); % time series of pp(<0.05Hz)
    for i = 1:nfreq % frequency from 0~5Hz
        ca_pp{ir}(:,:,i) = abs(tf2(i-1<freq2&freq2<i,:)); 
        ca_pp_ave{ir}(:,i) = squeeze(mean(ca_pp{ir}(:,:,i),1)); % timecourse pp(1~5Hz) individually
    end
end



% for ir = 1:length(scans)
%     h = figure;
%     % Q1: should 'ntimesout' be chosen to same samples as fmri?
%     % Q2: How does 'winsize' influence frequency resolution, indenpendent of NFFTs? 
%     [tf,freq,time] = timefreq(ca_demean{ir},fs_ca,'wletmethod','dftfilt3','winsize',win,'ntimesout',(len_ca/win*2),'padratio',256,'freqs',[0,10]);
%     time = time / 1000;
%    
%     imagesc(time,freq,abs(tf));
%     colormap jet;
%     set(gca,'YDir','normal');
%     xlabel('time(s)');
%     ylabel('Frequency(Hz)');
%     ylim([freq(1),freq(end)]);
%     set(gcf,'Position',[500 500 980 300 ]);
%     title(['ca 0.01-10Hz PSD trail ',num2str(scans(ir))]);
%     hcb = colorbar;
% 
%     caxis([0,6]);
%     ylabel(hcb,'Power(dB)');
% 
%     saveas(h,['ca 0.01-10Hz PSD trail ',num2str(scans(ir)),'.jpg']);
% %   
%     ca_pp_bsl{ir} = abs(tf(freq<0.05,:));% baseline assignment(<0.05Hz)
%     ca_pp_bls_ave{ir} = mean(ca_pp_bsl{ir},1); % time series of pp(<0.05Hz)
%     for i = 1:nfreq % frequency from 0~5Hz
%         ca_pp{ir}(:,:,i) = abs(tf(i-1<freq&freq<i,:)); 
%         ca_pp_ave{ir}(:,i) = squeeze(mean(ca_pp{ir}(:,:,i),1)); % timecourse pp(1~5Hz) individually
%     end
% 
% end
% Pwelch 

for ir = 1:length(scans)
% for ir = trails
    [pxx(:,ir),f] = pwelch(ca_demean{ir},win*2,[],[],fs_ca); 
%     plot(f(f>0.01&f<10),10*log10(pxx(f>0.01&f<10,ir)));
%     hold on;
end

% PSD plot [error bar]

h = figure;
pxx_min = min(10*log10(pxx(f>0.01&f<10,:)),[],2);
pxx_max = max(10*log10(pxx(f>0.01&f<10,:)),[],2);
k = fill([f(f>0.01&f<10);flipud(f(f>0.01&f<10))],[pxx_min;flipud(pxx_max)],[.9 .9 .9],'linestyle','none');
line(f(f>0.01&f<10),10*log10(mean(pxx(f>0.01&f<10,:),2)),'Color','black','LineWidth',2);
legend('all trails','mean of all trails');
xlabel('Frequency(Hz)');
ylabel('Power Density(dB/Hz)');
title('Averaged Calcium PSD');
box off;
saveas(h,'error bar PSD.jpg');
% saveas(h,['error bar PSD trail# ',num2str(ir),'.jpg']);
 for ir = 1:length(scans)
     h=figure;
     plot(f(f>0.01&f<10),10*log10(pxx(f>0.01&f<10,ir)),'LineWidth',2);
     xlim([0,10]);
     title(['Individual PSD trail# ',num2str(scans(ir))]);
     saveas(h,['Individual PSD trail# ',num2str(scans(ir)),'.jpg']);
 end
 close all
%% filtering calcium(0.01~0.1Hz)
% filtering ca(0.01~0.1Hz) timecourse
disp('ca filtering...');
bpf_low = 0.01;
bpf_high = 0.1;
[b,a] = butter(1,[bpf_low,bpf_high]/(fs_ca/2),'bandpass');
len = length(ca_demean{ir});
for ir = 1:length(scans)
freq_band_tmp = [fliplr(ca_demean{ir}(1:round(len/4))); ca_demean{ir}; ...
        fliplr(ca_demean{ir}(end-round(len/4)+1:end))];% extend signal for filtering
    ca_filt{ir} = filtfilt(b, a, freq_band_tmp);
    ca_filt{ir} = ca_filt{ir}(round(len/4)+1 : end-round(len/4));% trimmed  ca
end
% power profile 0.01~0.1Hz filtering
fs_filt_ca = 10;
[b1,a1] = butter(1,[bpf_low,bpf_high]/(fs_filt_ca/2),'bandpass');
len_filt = size(ca_pp_ave{ir},1);
clear freq_band_temp
% average of 1~5Hz
for ir = 1:length(scans)
ca_pp_all{ir} = mean(ca_pp_ave{ir}(:,2:end),2); % ave of 1~5Hz
end

% filter of 0~5Hz separately 
for ir = 1:length(scans)
    for i = 1:nfreq% 0~5 Hz frequency bands
        freq_band_tmp = [fliplr(ca_pp_ave{ir}(1:round(len_filt/4),i)); ca_pp_ave{ir}(:,i); ...
                fliplr(ca_pp_ave{ir}(end-round(len_filt/4)+1:end,i))];% extend signal for filtering
            ca_pb_filt{ir}(:,i) = filtfilt(b1, a1, freq_band_tmp);
            ca_pb_filt_new{ir}(:,i) = ca_pb_filt{ir}(round(len_filt/4)+1 : end-round(len_filt/4),i);% trimmed  ca
    end
end

% filter of 0~5Hz all
for ir = 1:length(scans)
    freq_band_tmp = [fliplr(ca_pp_all{ir}(1:round(len_filt/4),:)); ca_pp_all{ir}; ...
        fliplr(ca_pp_all{ir}(end-round(len_filt/4)+1:end,:))];% extend signal for filtering
    ca_pb_filt_all{ir} = filtfilt(b1, a1, freq_band_tmp);
    ca_pb_filt_all_new{ir} = ca_pb_filt_all{ir}(round(len_filt/4)+1 : end-round(len_filt/4),:);% trimmed  ca
end
clear freq_band_temp
% filter of <0.05Hz power
for ir = 1:length(scans)
    freq_band_tmp = [fliplr(ca_pp_bls_ave{ir}(:,1:round(len_filt/4)))'; ca_pp_bls_ave{ir}'; ...
        fliplr(ca_pp_bls_ave{ir}(:,end-round(len_filt/4)+1:end))'];% extend signal for filtering
    ca_pp_bls_filt{ir} = filtfilt(b1, a1, freq_band_tmp);
    ca_pp_bls_filt_new{ir} = ca_pp_bls_filt{ir}(round(len_filt/4)+1 : end-round(len_filt/4),:);% trimmed  ca
end
close all
%% Ca 1D plot

% ca time course(0.01~0.1Hz) 
for ir = 1:length(scans)
    h = figure;
    plot(t_ca,zscore(ca_filt{ir}));
    xlim([0,t_ca(end)]);
    xlabel('time(s)');
    ylabel('ca a.u.');
    set(gcf,'Position',[500 500 980 300 ]);
    title(['ca 0.01-0.1Hz trail ',num2str(ir)]);
%     saveas(h,['ca 0.01-0.1Hz timecourse trail ',num2str(ir),'.jpg'])
end
% power profile from PSD

% Ca freq-dependent power profile
for ir = 1:length(scans)
h = figure;
    for i = 1:5
        subplot(5,1,i);
        plot(1:length(time2),zscore(ca_pp_ave{ir}(:,i)));
%         xlim([time2(1),time2(end)]); [Modi.1]

        title(['Ca power profile ',num2str(i-1),'~',num2str(i), ' (Hz)']);
        mean_power(:,i) = mean(ca_pp_ave{ir}(:,i),1);
        txt = ['\mu(dB) = ', num2str(mean_power(:,i))];
        legend(txt,'FontSize',12,'TextColor','red','FontSize',10);
        legend('boxoff');
    end
%     saveas(h,['ca power profile trail ',num2str(ir),'.jpg']);
end
% Ca Power profile (0.01~0.1Hz)[separate]
for ir = 1:length(scans)
    h = figure;
    for i = 1:5
        subplot(5,1,i);
        plot(1:len_filt,zscore(ca_pb_filt_new{ir}(:,i)));
        xlim([0,640]);
        ylabel('Power(normalized)');
        title(['Filtered Ca power profile ',num2str(i-1),'~',num2str(i), ' (Hz)']);
    end
    saveas(h,['Filtered Ca power profile ',num2str(i-1),'~',num2str(i), ' (Hz).jpg']);
end
% Ca power profile (0.01~0.1Hz)[all averaged]
for ir = 1:length(scans)
    h = figure;
    plot(1:length(time2),zscore(ca_pb_filt_all_new{ir}));
    xlim([0,len_filt]);
    xlabel('time(s)');
    ylabel('Power(normalized)');
    title(['Filtered Ca power profile (0.01~10Hz) trail# ',num2str(ir)]);
    box off;
    set(gcf,'Position',[500 500 980 300 ]);
    saveas(h,['Filtered Ca power profile trail#',num2str(scans(ir)),'.jpg']);
end
close all
%% fmri filtering from 0.01-0.1Hz
[b2,a2] = butter(1, [0.01, 0.1]/((1/TR)/2), 'bandpass');
len_fmri = size(fmri{ir}{is},2);
for ir = 1:length(scans)
    for is = 1:nslice
        fmri_ave{ir}{is} = mean(fmri_tc{ir}{is},1);
        freq_band_tmp = [fliplr(fmri_ave{ir}{is}(1:round(len_fmri/4))'); fmri_ave{ir}{is}'; ...
            fliplr(fmri_ave{ir}{is}(end-round(len_fmri/4)+1:end)')];% extend signal for filtering
        fmri_filt{ir}{is} = filtfilt(b2, a2, freq_band_tmp);
        fmri_filt{ir}{is} = fmri_filt{ir}{is}(round(len_fmri/4)+1 : end-round(len_fmri/4));
    end
end


%
% voxel specific filter
for ir = 1:length(scans)
    for is = 1:nslice
        freq_band_tmp = [fliplr(fmri_trim{ir}{is}(:,1:round(len_fmri/4)))'; fmri_trim{ir}{is}'; ...
        fliplr(fmri_trim{ir}{is}(:,end-round(len_fmri/4)+1:end))'];% extend signal for filtering
        test222{ir}{is} = filtfilt(b2, a2, freq_band_tmp);
        test222{ir}{is} = test222{ir}{is}(round(len_fmri/4)+1 : end-round(len_fmri/4),:)';
    end
end
%

% filtered fmri(0.01~0.1Hz)
for ir = 1:length(scans)
    for is = 1:nslice
        h = figure;
        plot(t_fmri,zscore(fmri_filt{ir}{is}),'r');
        xlim([0,t_fmri(end)]);
        ylim([-3,3]);
        xlabel('time(s)');
        ylabel('fmri a.u.');
        set(gcf,'Position',[500 500 980 300 ]);
        title(['fmri 0.01-0.1Hz trail ',num2str(scans(ir)),' slice ',num2str(is)]);
        box off;
        saveas(h,['fmri 0.01-0.1Hz trail ',num2str(scans(ir)),' slice ',num2str(is),'.jpg']);
    end
end
close all
%% Ca & fMRI & Ca power profile
for is = 1:nslice
    for ir = 1:length(scans)
        h = figure;
        
        subplot(311);
        plot(1:length(time2),zscore(ca_pb_filt_all_new{ir}),'LineWidth',2,'Color',[0 0.4470 0.7410]);
        ylabel('Ca^2^+(1~5Hz)','fontweight','bold');
        xlim([0,6400]);
        box off;
        title(['Ultra-slow filtered signal trail# ',num2str(ir)]);
        
        subplot(312);
        plot(t_fmri,zscore(fmri_filt{ir}{is}),'LineWidth',2,'Color',[0.8500 0.3250 0.0980]);
        ylabel('fMRI (a.u.)','fontweight','bold');
        xlim([0,640]);
        box off;
        
        subplot(313);
        h1=plot(t_ca,ca_filt{ir},'LineWidth',2,'Color',[0.3010 0.7450 0.9330]);
        set(gcf,'Position',[500 500 980 300 ]);
        xlabel('time(s)');
        ylabel('Ca^2^+ time','fontweight','bold');
        xlim([0,640]);
        box off;
        saveas(h,['fmri_ca(1~5Hz) trail',num2str(ir),'slice',num2str(is),'.jpg']);
    end
    
end
close all
%% xcorr fmri slice
% cross slice fmri correlations
for ir = 1:length(scans)
    [fmri_corr12(ir,:),lag_0] = xcorr(zscore(fmri_filt{ir}{1}),zscore(fmri_filt{ir}{2}),100,'coeff');
    [fmri_corr13(ir,:),lag_0] = xcorr(zscore(fmri_filt{ir}{1}),zscore(fmri_filt{ir}{3}),100,'coeff');
    [fmri_corr23(ir,:),lag_0] = xcorr(zscore(fmri_filt{ir}{2}),zscore(fmri_filt{ir}{3}),100,'coeff');
end
% assign slice xcorr
xcor_fmri(:,:,1) = fmri_corr12;
xcor_fmri(:,:,2) = fmri_corr13;
xcor_fmri(:,:,3) = fmri_corr23;
% XCORR between 3 slices
for i = 1:3
err_min = min(xcor_fmri(:,:,i),2);
err_max = max(xcor_fmri(:,:,i),2);
h = figure;
fill([lag_0'*TR;flipud(lag_0'*TR)],[err_min';flipud(err_max')],[.9,.9,.9],'LineStyle','none');
line(lag_0*TR,mean(xcor_fmri(:,:,i),1),'Color','black','LineWidth',2);
[~,idx(:,i)] = max(mean(xcor_fmri(:,:,i),1),[],2);
xlabel('Lag time(sec)');
ylabel('correlation coefficient');
switch i
    case 1
        name_str = 'Cross-correlation Slice 1 & 2';
    case 2
        name_str = 'Cross-correlation Slice 1 & 3';
    case 3
        name_str = 'Cross-correlation Slice 2 & 3';   
        end
title(name_str);
xline((idx(:,i)-101)*TR,'LineWidth',1.5,'LineStyle','--');
box off;
legend('all trails','mean of all trails',['Peak Lag = ',num2str((idx(:,i)-101)*TR),' sec']);
saveas(h,[name_str,'.jpg']);
end
close all
% dowsample Ca
% ca_down{ir} = downsample(ca_pb_filt_new{ir},fs_ca*TR);
%% xcorr(fMRI & Ca time-course)
% if evoke == 1
    lg = 10;
    lg2 = 100;
% else
%     lg1 = 30;
%     lg2 = 300;
% end

for is = 1:nslice
      figure;
    for ir = 1:length(scans)
        ca_down(ir,:) = downsample(ca_filt{ir},fs_ca*TR);
        [xcf(ir,is,:),lag_1] = xcorr(zscore(ca_down(ir,:)),zscore(fmri_filt{ir}{is}),lg2,'coeff');% calcium leading fmri by how long
        plot(lag_1*TR,squeeze(xcf(ir,is,:)));
        hold on;
        xlabel('Lag time(sec)');
        ylabel('correlation coefficient');
        title(['cross-correlation(fmri&ca) slice',num2str(is)]);
        legend;
        [pks(ir,is,:),locs(ir,is,:)] = max(xcf(ir,is,:));
        locs(ir,is,:) = (locs(ir,is,:)-(lg2+1))*TR;
        plot(locs(ir,is,:),pks(ir,is,:),'*');
    end
end
%% average Ca^2^+ frequency power-band
for ir = 1:length(scans)
    ca_pb_filt_xc{ir} = mean(ca_pb_filt_new{ir},2);%1~5Hz [sperated]
end
%% xcorr (fMRI & Ca power-profile)


% xcorr (fmri & ca_power profile) [averaged] 
for is = 1:nslice
      h=figure;
      H = gca;
    for ir = 1:length(scans)
%         test_fmri(ir,is,:) = downsample(fmri_filt{ir}{is},1/TR);
        [xcf_test(ir,is,:),lag_2] = xcorr(zscore(ca_pb_filt_xc{ir}),zscore(squeeze(fmri_filt{ir}{is})),lg2,'coeff');% for 0-5Hz band
        lag = [-lg:lg/lg2:lg];
        h1= plot(lag,squeeze(xcf_test(ir,is,:)),'Color',[.65,.65,.65],'LineWidth',2);
% plot(lag_2,squeeze(xcf_test(ir,is,:)));
        hold on;
        xlabel('Lag time(sec)');
        ylabel('correlation coefficient');
        title(['Averaged Ca^2^+(1~5Hz) & fMRI corr slice',num2str(is)]);
        [pks_pp(ir,is,:),locs_pp(ir,is,:)] = max(xcf_test(ir,is,:));
        locs_pp(ir,is,:) = ((locs_pp(ir,is,:)-(lg2+1)))/lg;
        if locs_pp(ir,is,:)<0 || locs_pp(ir,is,:)==0 
            clr = 'r';
        else
            clr = 'b';    
        end
        plot(locs_pp(ir,is,:),pks_pp(ir,is,:),[clr,'*'],'MarkerSize',7);
        ylim([-1,1]);
        xticks([-10,-5,-3,-2,-1,0,1,2,3,5,10]);
        
        box off;
        set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%         legend('Positive lag','Negative lag(Ca leading)');
    end
    H.LineWidth = 1.5;
% saveas(h,['Negative correlated trails slice# ',num2str(is),'.jpg']);
saveas(h,['xcorr slice',num2str(is),'.jpg']);
end
% voxel-specific xcorr
for ir = 1:length(scans)
    for is = 1:nslice
     test333{ir}{is} = downsample(test222{ir}{is}',1/TR)';
    end
end

for ir = 1:length(scans)
    for is = 1:nslice
        for iv = 1:nvoxel
         [coeff_vo{ir}{is}(:,iv),lag_2] = xcorr(zscore(ca_pb_filt_xc{ir}),zscore(test222{ir}{is}(iv,:)'),lg2,'coeff');
        end
    end
end

% mean layer-specific lag time function handle
ave_lg = @(a,b) mean(a,2)-b;
for ir = 1:length(scans)
    for is = 1:nslice
        [~,lag_max_temp{ir}{is}] = max(coeff_vo{ir}{is},[],1);
        for il = 1:5
            lag_max{ir}{is}(:,il) = round(ave_lg(lag_max_temp{ir}{is}(:,voxel_idx==il),(lg2))/lg,1);
        end
        lag_max_tmp(ir,is,:) = lag_max{ir}{is};% pass cell to matrix
    end
end
% 1D xcorr [layer-specific]
Colour = colormap(hsv(5));
for ir = 1:length(scans)
    for is = 1:nslice
        h = figure;
        xlim([-10,10]);

        legd_names = {'L1','L2/3','L4','L5','L6'};
        for il = 1:5
            subplot(2,3,il);
            plot(lag,coeff_vo{ir}{is}(:,voxel_idx==il),'Color',Colour(il,:));
            ylabel([legd_names{il},' coeff']);
            xline(0,'k--');
            ylim([-1,1]);
            text(0,0.8,['\mu = ',num2str(lag_max{ir}{is}(:,il)),'s'],'Color',Colour(il,:));
            hold on;

        end
            
%         legend(legd_names');
        xlabel('lag time(s)');
        
        sgtitle(['Voxel-wise correlation trail# ',num2str(scans(ir)),' slice# ',num2str(is)]);
        saveas(h,['Voxel-wise correlation trail# ',num2str(scans(ir)),' slice# ',num2str(is),'.jpg']);
    end
end

% scatter layer-specific
for is = 1:nslice
    figure;
    for ir = 1:length(scans)
        plot(1:5,lag_max{ir}{is}/lg,'*');
        xticks([1,2,3,4,5]);
        hold on
    end
end
for is = 1:nslice
    for i = 1:5
        lg_ave(is,:,i) = mean(lag_max_tmp(:,is,i),1);
        lg_std(is,:,i) = std(lag_max_tmp(:,is,i),1);
    end
end
close all
%% error bar plot
h = figure;
barwitherr(squeeze(lg_std)',[1:5],squeeze(lg_ave)','LineWidth',2,'BarWidth',0.5);
legend({'Slice#1','Slice#2','Slice#3'});
set(gca,'XTickLabel',{'L1','L2/3','L4','L5','L6'});
% set(gca,'YDir','reverse');
colormap summer;
title('Maximum lag time (Ca&BOLD) in each layers');
box off;
grid on;
set(gca,'YLim',([-lg,lg]));
ylabel('Lag time(s)');
saveas(h,'Error bar.jpg');

%% 2D xcorr
for ir = 1:length(scans)
    h = figure;
    for is = 1:nslice
        test444{ir}(:,is) = mean(coeff_vo{ir}{is}(lag==0,:)',2);
    end
    imagesc(test444{ir});
    colormap jet;
    colorbar;
    title(['coefficient map at  lag trail# ',num2str(scans(ir))]);
    saveas(h,['coefficient map at  lag trail# ',num2str(scans(ir)),'.jpg']);
end   

for ir = 1:length(scans)
    h = figure;
    for is = 1:nslice
        subplot(1,3,is)
        test555{ir}(:,:,is) = coeff_vo{ir}{is};
        imagesc(lag,1:nvoxel,test555{ir}(:,:,is)');
    end
    sgtitle(['lag time map at  lag trail# ',num2str(scans(ir))]);
    saveas(h,['lag time map at  lag trail# ',num2str(scans(ir)),'.jpg']);
end   

% xcorr (fmri & ca_power profile) in [averaged] error map
for is = 1:nslice
    h = figure;
    xcf_max = max(squeeze(xcf_test(:,is,:)),[],1);
    xcf_min = min(squeeze(xcf_test(:,is,:)),[],1);
    fill([lag_2';flipud(lag_2')],[xcf_min';flipud(xcf_max')],[0.3010 0.7450 0.9330],'LineStyle','none');
    line(lag_2,squeeze(mean(xcf_test(:,is,:),1)),'Color',[0,0,1],'LineWidth',2);
    [~,idx_1(is,:)] = max(mean(xcf_test(:,is,:),1));
    xline(idx_1(is,:)-11,'LineWidth',1.5,'LineStyle','--','Color',[0 0.4470 0.7410]);
    xlabel('Lag time(s)');
    ylim([-1,1]);
    ylabel('Correlation coefficient');
    title(['Correlation of Ca^2^+(1~5Hz) power profile & fMRI slice# ',num2str(is)]);
    box off;
    legend('All trails','mean of all trails',['Averaged Peak Lag = ',num2str(idx_1(is,:)-11),' sec']);
    saveas(h,['Correlation of Ca(0.01~10Hz) power profile to fMRI slice# ',num2str(is),'.jpg']);
end

% fre-dependent xcorrelation (Ca & fMRI)

for i = 1:5
  h = figure;
    for ir = 1:length(scans)

        is = 3;

        [coff(ir,:,i),lag_3] = xcorr(zscore(ca_pb_filt_new{ir}(:,i)),zscore(squeeze(fmri_filt{ir}{is})),lg2,'coeff');
        h1 = plot(lag,coff(ir,:,i),'Color','k');
%         h1 = plot(lag_3,coff(ir,:,i));
        hold on;
        xlim([-lg,lg]);
        [pks_pp_new(ir,i,:),locs_pp_new(ir,i,:)] = max(coff(ir,:,i));
        locs_pp_new(ir,i,:) = (locs_pp_new(ir,i,:)-(lg2+1));
        if locs_pp_new(ir,i,:)<0 || locs_pp_new(ir,i,:)==0 
            clr = 'r';
        else
            clr = 'b';    
        end
        plot(locs_pp_new(ir,i,:)/lg,pks_pp_new(ir,i,:),[clr,'*'],'MarkerSize',7);
        ylim([-1,1]);
        box off;
        set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        legend('Positive lag','Negative lag(Ca leading)');
%         legend;
    end
    text(0.8,9,'slice 1');
    title(['Ca^2^+(',num2str(i-1),'~',num2str(i),'Hz) power band & fMRI correlation']);
    saveas(h,['Ca^2^+(',num2str(i-1),'~',num2str(i),'Hz) power band & fMRI correlation.jpg']);
end
close all
% trails = [5,6,7,8,9,10,12,14,15,16,17,18,19,20,21];%05232019
trails = [1,2,3,4,5,6];
clear  locs_pp_new
% freq dependent correlation 1D
% for ir = 1:length(scans)
% h = figure;
% for ir = 1:length(scans)
%     plot(1:5,locs_pp_new(ir,:,:),'*-','Color',[.9,.9,.9],'MarkerEdgeColor','r','LineWidth',2);
%     xlabel('Frequency(Hz)');
%     ylabel('lag time(s)');
%     hold on;
%     title('trails with similar correlation feature');
%     box off;
% end
% yline(0,'--','LineWidth',1);
% saveas(h,'corr_freq.jpg');

% trails = [2,3,4,5,7,9,10]; % 10062018[coupled]
% xcorr (2~3Hz & <0.05Hz)
h = figure;
for ir =1 : length(scans)
% for ir = trails
    [coeff_1(:,ir),lag] = xcorr(zscore(ca_pp_bls_filt_new{ir}),zscore(ca_pb_filt_all_new{ir}),lg2,'coeff');
    [pks_pb_new(ir,:),locs_pb(ir,:)] = max(coeff_1(:,ir));
    locs_pb(ir,:) = (locs_pb(ir,:)-(lg2+1));
    plot(lag,coeff_1(:,ir),'Color',[0.75 0.75 0.75]);
% plot(lag,coeff_1(:,ir));
    hold on;
    plot(locs_pb(ir,:),pks_pb_new(ir,:),'r*');
    xlabel('lag time(s)');
    ylabel('coeff');
    box off;
    title('xcorr between baseline & 1~5Hz Ca^2^+');
end
text(20,0.4,['N = ',num2str(ir)],'FontSize',14);
saveas(h,'xcorr between baseline & 1~5Hz Ca^2^+ .jpg');
