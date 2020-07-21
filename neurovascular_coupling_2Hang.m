%%%%%%%%%%%%%%%%%%%%%% USAGE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script aims to operate the following functions:
% - frequency dependent cross-correlation between calcium and fmri;(<0.05Hz & 2-3Hz Ca)
% - Statistical plot of BOLD along cortical depth;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fmri_path = 'D:\05232019\05232019_fmri_analyzed';
% fmri_path = 'D:\05302019_Hang\05302019_analyzed';
fmri_path = 'D:\03242019\fmri_analyzed03242019';
% fmri_path='D:\HANG\10062018rsesting-state\fmri_data';
% ca_path = 'D:\05232019\05232019_ca';
% ca_path = 'D:\05302019_Hang\05302019_ca';
ca_path = 'D:\03242019\Calcium_data03242019';
% ca_path = 'D:\HANG\10062018rsesting-state\calcium_data';
% mkdir('C:\Users\wangqi\Desktop\Neurovascular coupling\0324\figs');
savepath = 'C:\Users\wangqi\Desktop\Neurovascular coupling\0324\figs';
addpath('C:\Users\wangqi\Documents\Lab\Demo\Calcium');

TR = 0.1;% 0.1s
Spatial_res = 0.05;
prestim = 10;
fmri_duration = 640;

scans = [35 37 39 44 49 51 57 59 64 66 69];%03242019
% scans = [36,51,55,57,59,61,63,65,69,71];%10062018
% scans = [15 19 22 24 26 28 30 32 34 37 38 40 42 44 46 48 50 52 54 56 58 60 62 64 66 68 70 72 74 76 78 80];%05302019
% scans = [17 19 21 22 24 26 27 32 33 36 38 40 42 44 46 48 50 51 53 55 57];%05232019
fmri_dummy = ones(fmri_duration/TR,1);
chan_ca = 7;
nslice = 3;
cortical_depth = 2/Spatial_res;
offset = 1;

nfreq = 5;
%% load signals
% load calium
for ir = 1:length(scans)
    data_ca_orig = load([ca_path,'\scan_',num2str(scans(ir)),'.mat']);
    [data_match,fmri_dummy,beg,fin] = match_acq_fmri(data_ca_orig,fmri_dummy,TR,prestim);
    ca_match{ir} = -(data_match.channels{chan_ca}.data)';
end
fs_ca =  data_match.channels{chan_ca}.samples_per_second;
t_ca = [0:1/fs_ca:(length(ca_match{ir})-1)/fs_ca];
t_fmri = [0:TR:(length(fmri_dummy)-1)*TR];
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
% fmri trimming
% for is = 1:nslice
%     fmri_ccat{is} = cell2mat(fmri{is}');
%     mean_slice = mean(fmri_ccat{is}, 2);
%     mean_slice = mean_slice - min(mean_slice);
%     half_intens = max(mean_slice)/2; % half amplitude( thru out voxels)
%     i_half_ccat(is) = find(mean_slice>half_intens, 1);% half index(cortical voxel index)
% end
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


ir = 5;
is = 3;
[b1,a1] = butter(2,0.01/((1/TR)/2),'low');
freq_band_tmp = [fliplr(fmri_tc{ir}{is}(:,1:round(len_fmri/4))'); fmri_tc{ir}{is}'; ...
            fliplr(fmri_tc{ir}{is}(:,end-round(len_fmri/4)+1:end)')];% extend signal for filtering
        fmri_filt{ir}{is} = filtfilt(b1, a1, freq_band_tmp);
        fmri_filt{ir}{is} = fmri_filt{ir}{is}(round(len_fmri/4)+1 : end-round(len_fmri/4));


test = fmri_tc{ir}{is} - fmri_filt{ir}{is}'; % filtered off <0.01Hz
[b2,a2] = butter(1,[0.01,0.1]/((1/TR)/2),'bandpass');
freq_band_tmp = [fliplr(test(:,1:round(len_fmri/4))'); test'; ...
            fliplr(test(:,end-round(len_fmri/4)+1:end)')];% extend signal for filtering
        test_filt = filtfilt(b2, a2, freq_band_tmp);
        test_filt = test_filt(round(len_fmri/4)+1 : end-round(len_fmri/4));


test1 = downsample(test_filt,1/TR);
[pxx_test,lag] = xcorr(zscore(ca_pb_filt_new{5}(:,5)),zscore(test1),10,'coeff');
plot(lag,pxx_test);
xlabel('lagtime(s)');
ylim([-1,1]);

plot(1:640,test1,1:640,ca_pb_filt_new{5}(:,4))






% fmri time course each slice
% for ir = 1:length(scans)
for ir = 5
    for is = 1:nslice
        h = figure;
        plot(t_fmri,zscore(fmri_tc{ir}{is}),'r-');
        xlim([0,t_fmri(end)]);
        ylim([-3,3]);
        xlabel('time(s)');
        ylabel('fMRI a.u.');
        set(gcf,'Position',[500 500 980 300 ]);
        title(['fmri trail ',num2str(ir),' slice ',num2str(is)]);
        box off;
%         saveas(h,['fmri raw trail ',num2str(ir),' slice ',num2str(is),'.jpg']);
    end
end
% demean ca itself
for ir = 1:length(scans)
        ca_demean{ir} = ca_match{ir} - mean(ca_match{ir});
        h = figure;
%         plot(t_ca,zscore(ca_demean{ir}));
        plot(t_ca,ca_demean{ir});
        xlim([0,t_ca(end)]);
        xlabel('time(s)');
        ylabel('Ca^2^+ \DeltaF/F');
        set(gcf,'Position',[500 500 980 300 ]);
        title(['Ca^2^+ timecourse trail ',num2str(ir)]);
%         saveas(h,['Ca timecourse trail ',num2str(ir),'.jpg']);
end
%% spectral analysis(spectrogram & timefreq)
% ca spectrogram
disp('psd ca...')
win = 2*fs_ca;%2sec
NFFT = 2^18;
overlap = 0;
% for ir = 1:length(scans)
%     h = figure;
%     [S,F,T]= spectrogram(ca_demean{ir},blackman(win),overlap,NFFT,fs_ca,'yaxis');
%     imagesc(T,F(F<20),abs(S(F<20,:)));
%     set(gca,'YDir','normal');
%     xlabel('time(s)');
%     ylabel('Frequency(Hz)');
%     set(gcf,'Position',[500 500 980 300 ]);
% end
len_ca = length(ca_demean{ir});
for ir = 1:length(scans)
%     for ir = 2
    h = figure;
    % Q1: should 'ntimesout' be chosen to same samples as fmri?
    % Q2: How does 'winsize' influence frequency resolution, indenpendent of NFFTs? 
%     [tf,freq,time] = timefreq(ca_demean{ir},fs_ca,'winsize',win,'ntimesout',len_ca/win*2,'nfreqs',200,'freqs',[0.01,10]);
    [tf,freq,time] = timefreq(ca_demean{ir},fs_ca,'ffttaper','none','winsize',win,'ntimesout',len_ca/win*2,'padratio',256,'freqs',[0,10]);
    time = time / 1000;
   
    imagesc(time,freq,abs(tf));
    colormap jet;
    set(gca,'YDir','normal');
    xlabel('time(s)');
    ylabel('Frequency(Hz)');
    ylim([freq(1),freq(end)]);
    set(gcf,'Position',[500 500 980 300 ]);
    title(['ca 0.01-10Hz PSD trail ',num2str(ir)]);
    hcb = colorbar;
%     hold on;
%     plot(t_ca,zscore(ca_filt{ir})+6,'r-','LineWidth',1);% deploy filtered ca to PSD
%     
%     cmax = 0.98*max(max(abs(tf)));
    caxis([0,6]);
    ylabel(hcb,'Power(dB)');

    saveas(h,['ca 0.01-10Hz PSD trail ',num2str(ir),'.jpg']);
%   
    ca_pp_bsl{ir} = abs(tf(freq<0.05,:));% baseline assignment(<0.05Hz)
    ca_pp_bls_ave{ir} = mean(ca_pp_bsl{ir},1); % time series of pp(<0.05Hz)
    for i = 1:nfreq % frequency from 0~5Hz
        ca_pp{ir}(:,:,i) = abs(tf(i-1<freq&freq<i,:)); 
        ca_pp_ave{ir}(:,i) = squeeze(mean(ca_pp{ir}(:,:,i),1)); % timecourse pp(1~5Hz) individually
    end

end
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
     title(['Individual PSD trail# ',num2str(ir)]);
     saveas(h,['Individual PSD trail# ',num2str(ir),'.jpg']);
 end
%% filtering calcium(0.01~0.1Hz)
% filtering ca(0.01~0.1Hz) timecourse
disp('ca filtering...');
bpf_low = 0.01;
bpf_high = 0.1;
[b,a] = butter(1,[bpf_low,bpf_high]/(fs_ca/2),'bandpass');
len = length(ca_demean);
for ir = 1:length(scans)
freq_band_tmp = [fliplr(ca_demean{ir}(1:round(len/4))); ca_demean{ir}; ...
        fliplr(ca_demean{ir}(end-round(len/4)+1:end))];% extend signal for filtering
    ca_filt{ir} = filtfilt(b, a, freq_band_tmp);
    ca_filt{ir} = ca_filt{ir}(round(len/4)+1 : end-round(len/4));% trimmed  ca
end
% power profile 0.01~0.1Hz filtering
fs_filt_ca = 1;
[b1,a1] = butter(1,[bpf_low,bpf_high]/(fs_filt_ca/2),'bandpass');
len_filt = size(ca_pp_ave{ir},1);
clear freq_band_temp
% average of 1~5Hz
for ir = 1:length(scans)
ca_pp_all{ir} = mean(ca_pp_ave{ir}(:,2:end),2); % ave of 1~5Hz
end

% filter of 0~10Hz separately 
for ir = 1:length(scans)
    for i = 1:5% 0~10 Hz frequency bands
        freq_band_tmp = [fliplr(ca_pp_ave{ir}(1:round(len_filt/4),i)); ca_pp_ave{ir}(:,i); ...
                fliplr(ca_pp_ave{ir}(end-round(len_filt/4)+1:end,i))];% extend signal for filtering
            ca_pb_filt{ir}(:,i) = filtfilt(b1, a1, freq_band_tmp);
            ca_pb_filt_new{ir}(:,i) = ca_pb_filt{ir}(round(len_filt/4)+1 : end-round(len_filt/4),i);% trimmed  ca
    end
end

% filter of 0~10Hz all
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
        plot(1:length(time),zscore(ca_pp_ave{ir}(:,i)));
        xlim([time(1),time(end)]);
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
        plot(1:640,zscore(ca_pb_filt_new{ir}(:,i)));
        xlim([0,640]);
        ylabel('Power(normalized)');
        title(['Filtered Ca power profile ',num2str(i-1),'~',num2str(i), ' (Hz)']);
    end
    saveas(h,['Filtered Ca power profile ',num2str(i-1),'~',num2str(i), ' (Hz).jpg']);
end
% Ca power profile (0.01~0.1Hz)[all averaged]
for ir = 1:length(scans)
    h = figure;
    plot(1:length(time),zscore(ca_pb_filt_all_new{ir}));
    xlim([0,640]);
    xlabel('time(s)');
    ylabel('Power(normalized)');
    title(['Filtered Ca power profile (0.01~10Hz) trail# ',num2str(ir)]);
    box off;
    set(gcf,'Position',[500 500 980 300 ]);
    saveas(h,['Filtered Ca power profile trail#',num2str(ir),'.jpg']);
end

%% fmri filtering from 0.01-0.1Hz
[b2,a2] = butter(4, [0.01, 0.1]/((1/TR)/2), 'bandpass');
len_fmri = size(fmri{ir}{is},2);
for ir = 1:length(scans)
    for is = 1:nslice
        fmri_ave{ir}{is} = mean(fmri_trim{ir}{is},1);
        freq_band_tmp = [fliplr(fmri_ave{ir}{is}(1:round(len_fmri/4))'); fmri_ave{ir}{is}'; ...
            fliplr(fmri_ave{ir}{is}(end-round(len_fmri/4)+1:end)')];% extend signal for filtering
        fmri_filt{ir}{is} = filtfilt(b2, a2, freq_band_tmp);
        fmri_filt{ir}{is} = fmri_filt{ir}{is}(round(len_fmri/4)+1 : end-round(len_fmri/4));
    end
end
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
        title(['fmri 0.01-0.1Hz trail ',num2str(ir),' slice ',num2str(is)]);
        box off;
        saveas(h,['fmri 0.01-0.1Hz trail ',num2str(ir),' slice ',num2str(is),'.jpg']);
    end
end
%% Ca & fMRI & Ca power profile
for is = 1:nslice
    for ir = 1:length(scans)
        h = figure;
        
        subplot(311);
        plot(1:length(time),zscore(ca_pb_filt_all_new{ir}),'LineWidth',2,'Color',[0 0.4470 0.7410]);
        ylabel('Ca^2^+(1~5Hz)','fontweight','bold');
        xlim([0,640]);
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

%% xcorr fmri slice
% cross slice fmri correlations
for ir = 1:length(scans)
    [fmri_corr12(ir,:),lag] = xcorr(zscore(fmri_filt{ir}{1}),zscore(fmri_filt{ir}{2}),100,'coeff');
    [fmri_corr13(ir,:),lag] = xcorr(zscore(fmri_filt{ir}{1}),zscore(fmri_filt{ir}{3}),100,'coeff');
    [fmri_corr23(ir,:),lag] = xcorr(zscore(fmri_filt{ir}{2}),zscore(fmri_filt{ir}{3}),100,'coeff');
end
% assign slice xcorr
xcor_fmri(:,:,1) = fmri_corr12;
xcor_fmri(:,:,2) = fmri_corr13;
xcor_fmri(:,:,3) = fmri_corr23;
% XCORR between 3 slices
for i = 1:3
err_min = min(xcor_fmri(:,:,i));
err_max = max(xcor_fmri(:,:,i));
h = figure;
fill([lag'*TR;flipud(lag'*TR)],[err_min';flipud(err_max')],[.9,.9,.9],'LineStyle','none');
line(lag*TR,mean(xcor_fmri(:,:,i),1),'Color','black','LineWidth',2);
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
% dowsample Ca
% ca_down{ir} = downsample(ca_pb_filt_new{ir},fs_ca*TR);
%% xcorr(fMRI & Ca time-course)
for is = 1:nslice
      figure;
    for ir = 1:length(scans)
        ca_down(ir,:) = downsample(ca_filt{ir},fs_ca*TR);
        [xcf(ir,is,:),lag_1] = xcorr(zscore(ca_down(ir,:)),zscore(fmri_filt{ir}{is}),300,'coeff');% calcium leading fmri by how long
        plot(lag_1*TR,squeeze(xcf(ir,is,:)));
        hold on;
        xlabel('Lag time(sec)');
        ylabel('correlation coefficient');
        title(['cross-correlation(fmri&ca) slice',num2str(is)]);
        legend;
        [pks(ir,is,:),locs(ir,is,:)] = max(xcf(ir,is,:));
        locs(ir,is,:) = (locs(ir,is,:)-300)*TR;
        plot(locs(ir,is,:),pks(ir,is,:),'*');
    end
end
%% average Ca^2^+ frequency power-band
for ir = 1:length(scans)
    ca_pb_filt_xc{ir} = mean(ca_pb_filt_new{ir},2);
end
%% xcorr (fMRI & Ca power-profile)
% trails = [2,3,7,8,10,15,16,17,18,19];% negative correlated trails
trails = [2,3,5,7,8,9,12,14,15,16,17,18,23,24,25,26,27,30,31];
% xcorr (fmri & ca_power profile) [averaged] 
for is = 1:nslice
      h=figure;
    for ir = 1:length(scans)
% for ir = trails;

%     for ir = trails
        test_fmri(ir,is,:) = downsample(fmri_filt{ir}{is},1/TR);
        [xcf_test(ir,is,:),lag_2] = xcorr(zscore(ca_pb_filt_xc{ir}),zscore(squeeze(test_fmri(ir,is,:))),10,'coeff');% for 0-5Hz band
        h1= plot(lag_2,squeeze(xcf_test(ir,is,:)),'Color',[.9,.9,.9]);
% plot(lag_2,squeeze(xcf_test(ir,is,:)));
        hold on;
        xlabel('Lag time(sec)');
        ylabel('correlation coefficient');
        title(['Averaged Ca^2^+(1~5Hz) & fMRI corr slice',num2str(is)]);
        [pks_pp(ir,is,:),locs_pp(ir,is,:)] = max(xcf_test(ir,is,:));
        locs_pp(ir,is,:) = (locs_pp(ir,is,:)-11);
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
        legend('Positive lag','Negative lag(Ca leading)');
    end
% saveas(h,['Negative correlated trails slice# ',num2str(is),'.jpg']);
saveas(h,['xcorr slice',num2str(is),'.jpg']);
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
% for ir = 5
        
        fmri_dw{ir} = downsample(fmri_filt{ir}{1},1/TR);
        [coff(ir,:,i),lag_3] = xcorr(zscore(ca_pb_filt_new{ir}(:,i)),zscore(squeeze(fmri_dw{ir})),10,'coeff');
        h1 = plot(lag_3,coff(ir,:,i),'Color','k');
%         h1 = plot(lag_3,coff(ir,:,i));
        hold on;
        xlim([-10,10]);
        [pks_pp_new(ir,i,:),locs_pp_new(ir,i,:)] = max(coff(ir,:,i));
        locs_pp_new(ir,i,:) = (locs_pp_new(ir,i,:)-11);
        if locs_pp_new(ir,i,:)<0 || locs_pp_new(ir,i,:)==0 
            clr = 'r';
        else
            clr = 'b';    
        end
        plot(locs_pp_new(ir,i,:),pks_pp_new(ir,i,:),[clr,'*'],'MarkerSize',7);
        ylim([-1,1]);
        box off;
        set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        legend('Positive lag','Negative lag(Ca leading)');
%         legend;
end
    text(0.8,9,'slice 1');
    title(['Ca^2^+(',num2str(i-1),'~',num2str(i),'Hz) power band & fMRI correlation']);
%     saveas(h,['Ca^2^+(',num2str(i-1),'~',num2str(i),'Hz) power band & fMRI correlation.jpg']);
end

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
    [coeff_1(:,ir),lag] = xcorr(zscore(ca_pp_bls_filt_new{ir}),zscore(ca_pb_filt_all_new{ir}),30,'coeff');
    [pks_pb_new(ir,:),locs_pb(ir,:)] = max(coeff_1(:,ir));
    locs_pb(ir,:) = (locs_pb(ir,:)-31);
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
