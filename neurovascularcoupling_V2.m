

%%%%%%%%%%%%%%%%%%%%%% USAGE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script aims to operate the following functions:
% - Frequency dependent cross-correlation between calcium and fmri;(<0.05Hz & 2-3Hz Ca)
% - Statistical plot of BOLD along cortical depth;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PATH_1%%%%%%%%%%%%%%%%%%%%
% fmri_path = '/big_data/qi1/05232019/05232019_analyzed';
% fmri_path = '/big_data/qi1/11242019/11242019/11242019_fmri_raw/1.Wp1/';
% fmri_path = '/big_data/qi1/11212019/11212019/11212019_fMRI_raw/1.Wm1/';
% fmri_path = '/big_data/qi/10092018/10092018fmri/10092018';
% fmri_path = '/big_data/qi1/10062018/fmri_analyzed10062018';
% fmri_path = '/big_data/qi1/05302019/05302019_analyzed';
fmri_path = '/big_data/qi1/09122019/fmri_analyzed09122019';
% ca_path = '/big_data/qi1/05302019/05302019_ca';
% ca_path = '/big_data/qi1/11242019/11242019/11242019_3slice_ca';
% ca_path = '/big_data/qi1/05232019/05232019_ca';
% ca_path = '/big_data/qi1/11212019/11212019/11212019_3slice_calicum';
% ca_path = '/big_data/qi1/10062018/calcium10062018';
ca_path = '/big_data/qi1/09122019/09122019_ca';
mkdir('/home/qiwang/Documents/Cache_matlab/0912tk/figs');
savepath = '/home/qiwang/Documents/Cache_matlab/0912tk/figs';
addpath('/home/qiwang/SignalProcessing/Utilities');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PARAs%%%%%%%%%%%%%%%%%%%%%%%
nepoch = 32;
nY = 200;
TR = 0.1;% 0.1s
Spatial_res = 0.05; %0.05 per pixel/voxel
prestim = 10;%1 sec
fmri_duration = 640;
epoch_dur = fmri_duration / nepoch;
bsl_fmri = 2;%sec

rs = 0; % 1 | 0; rs | evoke
ave =1 ;
fmri_dummy = ones(fmri_duration/TR,1);
chan_ca = 7;
nslice = 3;
nvoxel = 2/Spatial_res;
cortical_depth = 2/Spatial_res;
offset = 100;
nfreq = 5;
voxel_idx = [1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5];
legd_names = {'L1','L2/3','L4','L5','L6'};


% Fonts Global
set(0,'DefaultAxesFontSize',15);
set(0,'DefaultAxesFontWeight','bold');
set(0,'DefaultAxesLineWidth',1.5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PATH_2%%%%%%%%%%%%%%%%%%%%
% scans = [35 37 39 44 49 51 57 59 64 66 69];%03242019
% scans = [36,40,42,55,57,61,63,65,69];%10062018rs_GOOD
% scans = [30,36,38,40,42,43,45,47,49,51,53,55,57,59,61,63,65,67,69,71,73];%10062018_RS
% scans = [24,29,31,35,37,39,41,44,46,48,50,52,54,56,58,60,62,64,66,68,70,72,74];%10062018TK
% scans = [13,15,17,19,21,23,25,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64,66,68,70,72,74,76,78];%09122019_RS
scans = [12,14,16,18,20,22,24,27,29,31,33,35,37,39,41,43,45,47,49,51,53,55,57,59,61,63,65,67,69,71,73,75,77,79];%09122019_TK
% scans = [26,28,30,34,44,50,54,56,58,74,76];%09122019RS_GOOD
% scans = [12,14,15,16,18,20,22,26,39,41,43,45,47,49,51,53,55,57,59,61,63,65,67,69];%1121_evoked
% addpath('/big_data/qi1/11242019/11242019/11242019_3slice_ca/PSD_RS/psd_mat')
% addpath('/big_data/qi/05232019/05232019_ca/PSD_RS/psd_data')
% addpath('/big_data/qi1/11212019/11212019/11212019_3slice_calicum/PSD_RS/psd_mat')
% addpath('/big_data/qi1/05232019_psd/psd_tk');
% addpath('/big_data/qi1/05232019_psd')
% addpath('/big_data/qi1/05302019/psd_tk');
% addpath('/big_data/qi1/05302019/psd_rs');
% addpath('/big_data/qi1/10062018/calcium10062018/psd_RS/psd_mat');
% addpath('/big_data/qi1/10062018/calcium10062018/psd_TK/psd_mat');
addpath([ca_path,'/TK_ca_psd']);%0912
% addpath([ca_path,'/RS_ca_psd']);%0912rs


% scans = [13,17,19,21,23,25,27,29,40,42,44,46,48,50,52,54,56,58,60,62,64,66,68,70];%1121RS
% scans = [42,44,48,50,56,60,62,66,70];%1121rs_good
% scans = [11,21,23,25,28,41,43,45,47,49,51,53,55,57,59,61,63,64,66];%11242019_rs
% scans =  [55,57,59,61,64];%1124rs_good
% scans = [10,22,24,26,29,40,42,44,46,48,50,52,54,56,58,60,62,67];%1124_evoked
% scans = [16 17 20 23 25 27 29 31 33 35 36 ];%05302019tk_P
% scans = [39 41 43 45 47 49 51 53 55 57 59 61 63 65 67 69 71 73 75 77 79 81];%05302019tk_N
% scans = [39 41 43 45 47 49 51 53 55 57 59 61 63 65 67 69 71 73 75 77 79 81];%05302019tk
% scans = [15 19 22 24 26 28 30 32 34 37 38 40 42 44 46 48 50 52 54 56 58 60 62 64 66 68 70 72 74 76 78 80];%05302019_RS
% scans = [44,46,48,50,52,54,60,62,64,66,68,72,74,76];%0530rs_good
% scans = [17 19 21 22 24 26 27 32 33 36 38 40 42 44 46 48 50 51 53 55 57];%05232019RS
% scans = [15 18 20 23 25 28 31 34 35 37 39 41 43 45 52 54 56 58 59];%05232019TK
% scans = [26,33,40,42,44,48,50,55,57];%0523_rs_GOOD

 % good trials

g_scans =scans;%1124
 

%% load signal

% load calium
for ir = 1:length(scans)
    data_ca_orig = load([ca_path,'/scan',num2str(scans(ir)),'.mat']);
    [data_match,fmri_dummy,beg,fin] = match_acq_fmri(data_ca_orig,fmri_dummy,TR,prestim);
    ca_match{ir} = -(data_match.channels{chan_ca}.data)';
    fs_ca =  data_match.channels{chan_ca}.samples_per_second;
%     baseline_ca{ir} = mean(-data_ca_orig.channels{7}.data(:,1:1*fs_ca),2); 
end
% set time scales

t_ca = [0:1/fs_ca:(length(ca_match{ir})-1)/fs_ca];
t_fmri = [0:TR:(length(fmri_dummy)-1)*TR];
t_ep = [0:TR:(nY-1)*TR];
% load fmri
for ir = 1:length(scans)
    for is = 1:nslice
        tmp = load([fmri_path, '/', num2str(scans(ir)), '/Results_Slice', num2str(is), ...
                '/line_scanning_data_s', num2str(is), '.mat']);
        
        fmri_load = abs(tmp.total_cortical_depth_map);
%          for 1124 only
if strcmp(fmri_path, '/big_data/qi1/11242019/11242019/11242019_fmri_raw/1.Wp1/')
        fmri_raw(ir,4-is,:,:) = fmri_load;  
else
        fmri_raw(ir,is,:,:) = fmri_load; % raw fmri_input   
end
    end
end

if ave ==1
    trials = 1;
    fmri_raw = mean(fmri_raw,1); % ave raw
else
    trials = length(scans);
end
%% reshaping pre-process
cd (savepath);
% [DEFAULT SETTING]
% for ir = 1:trials
%     for is = 1:nslice
%         mean_slice = mean(fmri{ir}{is}, 2);
%         mean_slice = mean_slice - min(mean_slice);% mean of one scan
%         
%         half_intens = max(mean_slice)/2;
%         i_half(is,ir) = find(mean_slice>half_intens, 1);% the first voxel with value larger than trail average
%         
% %         shift the cortex lower by several voxels
% %         i_half(is,ir) = i_half(is,ir)+offset;
%     end
% end

% % % ONLY VALID FOR 1006RS:
% i_half = zeros(nslice,trials);
% 
% i_half(1,1:13) = 50;
% i_half(2,1:13) = 58;
% i_half(3,1:13) = 60;
% i_half(1,14:end) = 63;
% i_half(2,14:end) = 62;
% i_half(3,14:end) = 66;
% ONLY VALID FOR 0523RS:
% i_half = zeros(nslice,trials);
% 
% i_half(1,:) = 39;
% i_half(2,:) = 39;
% i_half(3,:) = 41;
% ONLY VALID FOR 0530RS:
% i_half = zeros(nslice,trials);
% 
% i_half(1,:) = 29;
% i_half(2,:) = 26;
% i_half(3,:) = 28;

% % ONLY VALID FOR 1121TK:
% i_half = zeros(nslice,trials);
% i_half(1,:) = 35;
% i_half(2,:) = 35;
% i_half(3,:) = 35;
% 
% i_half(1,1) = 40;
% i_half(1,2:8) = 55;
% i_half(1,9:end)= 35;
% i_half(2,1) = 40;
% i_half(2,2:8) = 55;
% i_half(2,9:end)= 35;
% i_half(3,1) = 40;
% i_half(3,2:8) = 55;
% i_half(3,9:end)= 36;
% trim 2mm of cortex

% % ONLY VALID FOR 1121RS:

% i_half = zeros(nslice,trials);
% 
% i_half(1,:) = 34;
% i_half(2,:) = 39;
% i_half(3,:) = 35;
% trim 2mm of cortex


% ONLY VALID FOR 1124TK:
% i_half = zeros(nslice,trials);
% 
% i_half(1,:) = 50;
% i_half(2,:) = 33;
% i_half(3,:) = 45;

% trim 2mm of cortex

% ONLY VALID FOR 0912RS:
i_half = zeros(nslice,trials);

i_half(1,:) = 43;
i_half(2,:) = 34;
i_half(3,:) = 34;
%% Trim fMRI
for ir = 1:trials
    for is = 1:nslice
        fmri_trim{ir}{is} = squeeze(fmri_raw(ir,is,i_half(is,ir)+1:i_half(is,ir)+cortical_depth,:));
        fmri{ir}{is} = squeeze(fmri_raw(ir,is,:,:));
    end
end

% tSNR & cortical map

for ir  = 1:trials
    for is  = 1:nslice
        h = figure;
        subplot(121);
        imagesc(t_fmri,1:size(fmri{ir}{is},1),fmri{ir}{is});
        hold on;
        colormap gray;
        xlabel('time(s)');
        ylabel('cortical depth');
        yline(i_half(is,ir),'b--','LineWidth',2);
        yline(i_half(is,ir)+cortical_depth,'b--','LineWidth',2);
        hold on;
        subplot(122);
        plot(mean(fmri{ir}{is},2));
        xline(i_half(is,ir),'b--','LineWidth',2);
        xline(i_half(is,ir)+cortical_depth,'b--','LineWidth',2);
        xlim([0,128]);
        sgtitle(['SNR trail ',num2str(scans(ir)),' slice ',num2str(is)]);
        set(gcf,'Position',[500 500 980 300 ]);
        saveas(h,['SNR trail ',num2str(scans(ir)),' slice ',num2str(is),'.jpg']);
    end
end
%assign ca & fmri
for ir = 1:trials
    for is = 1:nslice
%         fmri_tc{ir}{is} = detrend(mean(fmri_trim{ir}{is},1));
        fmri_tc{ir}{is} = mean(fmri_trim{ir}{is},1);
    end
end
% close all
%% Percentage change of signal

% Ca ave-trials
if rs == 0
for ir = 1:length(scans)
    for i = 1:nepoch
        prestim_ca{ir}(:,i) = mean(ca_match{ir}(1+(i-1)*fs_ca*epoch_dur:prestim*(fs_ca*TR)+(i-1)*fs_ca*epoch_dur));%prestim ca baseline
        ca_ep_tmp(:,i) = (ca_match{ir}(1+(i-1)*fs_ca*epoch_dur:i*fs_ca*epoch_dur) - prestim_ca{ir}(:,i))./prestim_ca{ir}(:,i);% percentage change
        ca_max_tmp(:,i) = max(ca_ep_tmp(:,i),[],1);
    end
    per_ca{ir} = reshape(ca_ep_tmp,1,[]);
    ca_ep{ir} = mean(ca_ep_tmp,2);
    ca_max{ir} = ca_max_tmp;
end
ca_ave_ep = mean(cell2mat(ca_ep),2);
end

% fMRI
% norm to 100
perc = @(x) (x*100)./mean(x,2);
for ir = 1:trials
    for is = 1:nslice
        fmri_norm{ir}{is} = perc(fmri_trim{ir}{is});% raw data norm
        fmri_de{ir}{is} = (detrend(fmri_norm{ir}{is}',2))';%6400*40
        per_fmri{ir}{is} = ((perc(fmri_de{ir}{is}+10) - 100)/100);%6400*40
    end
end
% epoch-wise averaging
if rs == 0
    for ir = 1:trials
        for is = 1:nslice
            for i = 1:nepoch
                prestim_fmri{ir}{is}(:,i) = mean(fmri_norm{ir}{is}(:,1+(i-1)*nY:prestim+(i-1)*nY),2);%prestim per epoch
                fmri_tmp(:,:,i) = (fmri_norm{ir}{is}(:,1+(i-1)*nY:i*nY) - prestim_fmri{ir}{is}(:,i)) ./prestim_fmri{ir}{is}(:,i); %each epoch
                fmri_max_tmp(:,i) = max(fmri_tmp(:,:,i),[],2);
            end
            per_fmri{ir}{is} = reshape(fmri_tmp,nvoxel,nY*nepoch);
            fmri_ep{ir}{is} = mean(fmri_tmp,3);
            fmri_max{ir}{is} = fmri_max_tmp;% peak of each epoch
        end
    end
end
%% Onset time & peak amplitude [fMRI][NOT distinguishing]
if ave == 1
for ir = 1:trials
    h=figure;
 for is = 1:nslice     
        for iv = 1:nvoxel
            prestim_std{ir}{is}(iv) = std(100*fmri_ep{ir}{is}(iv,1:prestim));
            [~,tt(iv,:)] = find(100*fmri_ep{ir}{is}(iv,prestim+1:end)>(1*prestim_std{ir}{is}(iv)),1);
%             figure;
%              plot(t_ep,100*fmri_ep{ir}{is}(iv,:));
%              hold on;
%              yline(2*prestim_std{ir}{is}(iv),'Color','red');
%              xline((prestim+tt(iv))*TR);
%              yline(100*fmri_ep{ir}{is}(iv,prestim+tt(iv)),'LineStyle','--');
        end
        subplot(3,1,is)
        [p1,~,mu1] = polyfit(1:nvoxel,tt',12);
        FF = polyval(p1,1:nvoxel,[],mu1);
        plot(1:nvoxel,FF*TR,'k-','LineWidth',2);
        ylim([0,2]);
        yticks([0,0.5,1,1.5,2]);
        if is ==3
        xlabel('cortical depth');
        elseif is ==2
        ylabel('On-set time (s)');        
        end
        onset{ir}{is}= tt';

 end
 saveas(h,['On set time.jpg']);

end

for ir = 1:trials
    h = figure;
    for is = 1:nslice
        subplot(3,1,is)
        [p,~,mu] = polyfit(1:nvoxel,prestim_std{ir}{is},12);% std 12-th polyfit
        F = polyval(p,1:nvoxel,[],mu);% fitted
        plot(1:nvoxel,F,'k-','LineWidth',2);
        ylim([0,0.5]);
        yticks([0,0.1,0.2,0.3,0.4,0.5])
        if is ==3
        xlabel('cortical depth');
        elseif is ==2
        ylabel('standard deviation');        
        end
    end
     saveas(h,['std of baseline.jpg']);
end
end

if ave ==1
h = figure;
for ir = 1:trials
    for is = 1:nslice
        plot(100*mean(fmri_max{ir}{is},1),100*ca_max{ir},'.','Color',[0.1+(is-1)*0.3,0.1+(is-1)*0.3,0.1+(is-1)*0.3],'DisplayName',['Slice#',num2str(is)],'MarkerSize',10); 
    %     xlim([0,10]);
        xlabel('BOLD(%)');
    %     ylim([30,80]);
        ylabel('Ca^2^+(%)');
        hold on;   
%         legend ;
    end
end
saveas(h,['amplitude correlation.jpg']);
end
%% Raw signals 1-D plot
cd (savepath)
for ir = 1:trials
    h = figure;
    subplot(4,1,1);
    plot(t_ca/60,100*per_ca{ir},'Color',[107,72,255]/255);
    ylabel('Ca^2^+ (\Delta F / F)');
    xlim([0,t_ca(end)/60]);
    box off;
    ax = gca;
    ax.XAxis.Visible = 'off';
    
    for is = 1:nslice
        subplot(4,1,is+1)
        plot(t_fmri/60,100*mean(per_fmri{ir}{is},1),'k-');
        if is ~= 3
            ax = gca;
            ax.XAxis.Visible = 'off';
        end
        xlim([0,t_fmri(end)/60]);
        ylim([-1,6]);
        xlabel('time(min)');
        ytickformat('%g%%');
%         ylim([-1,5]);
        ylabel({'BOLD (%)';['Slice #',num2str(is)]});
        set(gcf,'Position',[1070 500 1770 840 ]);
        box off;
    end
    saveas(h,['Raw timecourse trail ',num2str(scans(ir)),'.jpg']);
end

% fmri time course each slice

% 1D epoch fmri response
clr_ep = flip(colormap(jet(nvoxel)),1);
for ir = 1:trials
    h = figure;
    for is = 1:nslice
        subplot(1,3,is)
        for iv = nvoxel:-1:1
            plot(t_ep,100*fmri_ep{ir}{is}(iv,:),'Color',clr_ep(iv,:));
            hold on;
            xlabel('Times(s)');
            ytickformat('%g%%');
            ylabel('BOLD change');
        end
        if is ~= 1
        ax = gca;
        ax.YAxis.Visible = 'off';
        end
        xticks([0 10 20])
        ylim([-10,10]);
        box off;
        h.Position = [1040,1520,1670,300];
        if ave == 1
            title(['Slice#',num2str(is)]);
        else
            title(['Epoch BOLD trial#',num2str(scans(ir)),' slice#',num2str(is)]);
        end
    end
    colormap(clr_ep);
%     colorbar('Ticks',[0,1],'Direction','reverse','TickLabels',{'Surface','Bottom'},'Location','eastoutside');
    saveas(h,['Epoch BOLD trial#',num2str(scans(ir)),'.jpg']);   
end
% 2D epoch fmri
for ir = 1:trials
    for is = 1:nslice
        h = figure;
        imagesc(t_ep,[],100*fmri_ep{ir}{is});
        colormap jet;
        colorbar;
    %     caxis([-1,5]);
    for il = [4,12,17,28,40]
        yline(il,'LineStyle','--','Color',[1,1,0],'LineWidth',2.5);
        text(0,il-2,legd_names{voxel_idx(il-1)},'Color','white','FontSize',14,'FontWeight','bold');
    end
        xlabel('time(s)');
        ylabel('Cortical depth(mm)');
        yticks([0,10,20,30,40]+0.5);% 0.5 for imagesc() offset from ori tick
        yticklabels({'0','0.5','1','1.5','2.0'});
        saveas(h,['epoch 2d raw trials#',num2str(scans(ir)),' S#',num2str(is),'.jpg']);
    end
end

if ave == 0
    % CNR computation
    for ir = 1:trials
        for is = 1:nslice
            ma = max(fmri_ep{ir}{is},[],2);
            stdvv = std(fmri_ep{ir}{is}(:,1:prestim/TR),[],2);
            CNR(ir,is,:) = ma ./ stdvv;
        end
    end
    cnr_m = mean(CNR,1);
    cnr_std = std(CNR,1);
    for is = 1:nslice
        h = figure;
        plot(1:nvoxel,squeeze(cnr_m(:,is,:)),'Color','r','LineWidth',2.5);
        hold on;
        errorbar(1:nvoxel,squeeze(cnr_m(:,is,:)),squeeze(cnr_std(:,is,:)),'Color',[.7,.7,.7],'LineWidth',1);
        title(['CNR of slice#',num2str(is)]);
        xlabel('Cortical depth');
        ylabel('CNR');
        saveas(h,['CNR of slice#',num2str(is),'.jpg']);
    end
end
% ca 1D
% for ir = 1:trials
%         h = figure;
%         plot(t_ca,per_ca{ir}*100);
%         xlim([0,t_ca(end)]);
%         xlabel('time(s)');
%         ylabel('Ca^2^+ \Delta F/F');
%         ytickformat('%g%%');
%         set(gcf,'Position',[500 500 980 300 ]);
%         title(['Ca^2^+ timecourse trail ',num2str(scans(ir))]);
%         saveas(h,['Ca timecourse trail ',num2str(scans(ir)),'.jpg']);
% end
t_ca_ep = [0:1/fs_ca:(length(ca_ep{ir})-1)/fs_ca];
for ir = 1:trials
        h = figure;
        plot(t_ca_ep,ca_ep{ir}*100);
%         xlim([0,t_ca_ep(end)]);
        xlabel('time(s)');
        ylabel('Ca^2^+ \Delta F/F');
        ytickformat('%g%%');
%         set(gcf,'Position',[500 500 980 300 ]);
        title(['Ca^2^+ epoch trail ',num2str(scans(ir))]);
        saveas(h,['Ca epoch trail ',num2str(scans(ir)),'.jpg']);
end

h = figure;
plot(t_ca_ep,ca_ave_ep*100,'Color',[0,0,1],'LineWidth',1.5);
ytickformat('%g%%');
ylabel('\Delta F / F');
xlabel('time(s)');
% ylim([-2,20]);
xticks([0,5,10,15,20]);
% ax = gca;
% ax.YAxis.Visible = 'off'
box off;
% title('Ca^2^+ epoch\_ave ');
saveas(h,'Ca epoch_ave.jpg');
close all

% demean ca itself
for ir = 1:trials
ca_demean{ir} =  ca_match{ir} - mean(ca_match{ir});
end
%% spectral analysis(spectrogram & timefreq)
% ca spectrogram
disp('psd ca...')
win = 2*fs_ca;%2sec
NFFT = 2^18;
overlap = 0;
len_ca = length(ca_demean{ir});


% PSD producing
% cd([ca_path,'/psd_TK/psd_mat'])
% for ir = 1:trials
%     for c = 1:4
%         [tf1(:,:,c),freq1,time1(:,c)] = timefreq(ca_demean{ir}((len_ca/4)*(c-1)+1:c*(len_ca/4),:),fs_ca,'wletmethod','dftfilt3','winsize',win,'ntimesout',6400/4,'padratio',256,'freqs',[0,10]);
%     end
%     save(['tf',num2str(scans(ir)),'.mat'], 'tf1');
%     clear tf1
% end

for ir = 1:trials
    tf1{ir} = load(['tf',num2str(scans(ir)),'.mat']);
end
cd(savepath)
for ir = 1:trials
    h = figure;
tf2{ir} = reshape(tf1{ir}.tf1,size(tf1{ir}.tf1,1),6400);
freq2 = 0:10/(size(tf1{ir}.tf1,1)-1):(size(tf1{ir}.tf1,1)-1)*10/(size(tf1{ir}.tf1,1)-1); % set freq ticks array 
time2 = 0:.1:(6400-1)*.1;% set time ticks array
imagesc(time2,freq2,abs(tf2{ir}));
colormap jet;
caxis([0,6]);
set(gca,'YDir','normal');
xlabel('time(s)');
ylabel('Frequency(Hz)');
ylim([0,freq2(end)]);
set(gcf,'Position',[500 500 980 300 ]);
saveas(h,['ca 0.01-10Hz PSD trail ',num2str(scans(ir)),'.jpg']);
end
for ir = 1:trials
ca_pp_bsl{ir} = abs(tf2{ir}(freq2<0.05,:));% baseline assignment(<0.05Hz)
    ca_pp_bls_ave{ir} = mean(ca_pp_bsl{ir},1); % time series of pp(<0.05Hz)
    for i = 1:nfreq % frequency from 0~5Hz
        ca_pp{ir}(:,:,i) = abs(tf2{ir}(freq2==i-1|i-1<freq2&freq2<i,:)); 
        ca_pp_ave{ir}(:,i) = squeeze(mean(ca_pp{ir}(:,:,i),1)); % timecourse pp(1~5Hz) individually
    end
end



 

for ir = 1:trials
% for ir = trails
    [pxx(:,ir),f] = pwelch(ca_demean{ir},win*2,[],[],fs_ca); 
    [pxx_low(:,ir),f_low] = pwelch(ca_demean{ir},length(pwelch(ca_demean{ir})),[],[],fs_ca);
%     plot(f(f>0.01&f<10),10*log10(pxx(f>0.01&f<10,ir)));
%     hold on;
end

% PSD plot [error bar]
 %% arteriole spectra

% individual trail ca_psd & low oscillation
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
 for ir = 1:trials
     h=figure;
     subplot(311)
     plot(f(f>0.01&f<10),10*log10(pxx(f>0.01&f<10,ir)),'LineWidth',2);
     xlim([0,10]);
     ylabel('Power(10*log10 scale)');
     title(['Individual PSD trail# ',num2str(scans(ir))]);
     subplot(312)
     plot(f_low(f_low<0.1),pxx_low(f_low<0.1,ir),'LineWidth',2);
     ax = gca();
     ax.XScale = 'log';
     xticks([0,0.01,0.02,0.04,0.08,0.1]);
     xlim([f_low(1),0.1]);
     ylabel('Power(dB)');

%      xline(100*ls(ir,:),'--','LineWidth',3);
     subplot(313)
     plot(f_low(f_low<0.1),10*log10(pxx_low(f_low<0.1,ir)),'LineWidth',2);
     ax = gca();
     ax.XScale = 'log';
     xticks([0,0.01,0.02,0.04,0.08,0.16,0.3,0.5]);
     xlim([f_low(1),0.1]);
     xlabel('Frequency(Hz)');     
     ylabel('Power(10*log10 scale)');
%      xline(100*ls(ir,:),'--','LineWidth',3);
     
     saveas(h,['Individual PSD trail# ',num2str(scans(ir)),'.jpg']);
 end
 close all
 % fmri color palette
 fcp = [1,0,0;0,1,0;0,0,1];
 h = figure;
 for ir = 1:trials
     if any(g_scans==scans(ir))
%       h = figure;
    for is = 1:nslice
        
        subplot(3,1,is)
        [pxxx(:,ir,is),ff] = pwelch(mean(per_fmri{ir}{is},1),length(fmri_dummy),[],[],1/TR);
        plot(ff(ff<.12),log10(pxxx(ff<.12,ir,is)),'LineWidth',1,'Color',fcp(is,:) );
%         plot(diff(log10(pxxx(ff<.12,ir,is)))/0.001,'LineWidth',1,'Color',fcp(is,:) );
        hold on;
        box off;
        ax = gca();
%         ax.XScale = 'log';
%         xticks([0,0.005,0.01,0.02,0.05,0.1,0.3,0.5]);
        xlabel('Freq(Hz)');
        ylabel('Power(dB)');
%         xlim([ff(1),0.12]);
%         title(['Trial#',num2str(scans(ir)),' Slice#',num2str(is)]);
title(['slice#',num2str(is)]);
    end
    legend ({'slice#1','slice#2','slice#3'});
%     saveas(h,['fmri_psd trial#',num2str(scans(ir)),'.jpg']); 
     end
    
 end
 title('coupled trials PSD');
figure; for ir = 1:trials
     if any(g_scans==scans(ir))
%       h = figure;
    for is = 1:nslice
        
        subplot(3,1,is)
        [pxxx(:,ir,is),ff] = pwelch(mean(per_fmri{ir}{is},1),length(fmri_dummy),[],[],1/TR);
        plot(ff(ff<.12),log10(pxxx(ff<.12,ir,is)),'LineWidth',1,'Color',fcp(is,:) );
        hold on;
        box off;
        ax = gca();
%         ax.XScale = 'log';
%         xticks([0,0.005,0.01,0.02,0.05,0.1,0.3,0.5]);
        xlabel('Freq(Hz)');
        ylabel('Power(dB)');
%         xlim([ff(1),0.1]);
%         title(['Trial#',num2str(scans(ir)),' Slice#',num2str(is)]);
title(['slice#',num2str(is)]);
    end
    legend ({'slice#1','slice#2','slice#3'});
%     saveas(h,['fmri_psd trial#',num2str(scans(ir)),'.jpg']); 
     end
    
 end
%   for ir = 1:trials
%      if ~any(g_scans==scans(ir))
% %       h = figure;
%     for is = 1:nslice
%         
%         subplot(3,1,is)
%         [pxxx(:,ir,is),ff] = pwelch(per_fmri{ir}{is},length(fmri_dummy)/8,[],[],1/TR);
%         plot(ff(ff<.5),log10(pxxx(ff<.5,ir,is)),'LineWidth',1,'Color',fcp(is,:) );
%         hold on;
%         box off;
%         ax = gca();
% %         ax.XScale = 'log';
% %         xticks([0,0.005,0.01,0.02,0.05,0.1]);
%         xlabel('Freq(Hz)');
%         ylabel('Power(dB)');
% %         xlim([ff(1),0.1]);
%         title(['Trial#',num2str(scans(ir)),' Slice#',num2str(is)]);
%     end
%     legend ({'slice#1','slice#2','slice#3'});
% %     saveas(h,['fmri_psd trial#',num2str(scans(ir)),'.jpg']); 
%      end
% 
%  end
%% filtering calcium(0.01~0.1Hz)
% filtering ca(0.01~0.1Hz) timecourse
disp('ca filtering...');
bpf_low = 0.01;
bpf_high = 0.1;
[b,a] = butter(1,[bpf_low,bpf_high]/(fs_ca/2),'bandpass');
len = length(ca_demean{ir});
for ir = 1:trials
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
for ir = 1:trials
ca_pp_all{ir} = mean(ca_pp_ave{ir}(:,2:end),2); % ave of 1~5Hz
end

% filter of 0~5Hz separately 
for ir = 1:trials
    for i = 1:nfreq% 0~5 Hz frequency bands
        freq_band_tmp = [fliplr(ca_pp_ave{ir}(1:round(len_filt/4),i)); ca_pp_ave{ir}(:,i); ...
                fliplr(ca_pp_ave{ir}(end-round(len_filt/4)+1:end,i))];% extend signal for filtering
            ca_pb_filt{ir}(:,i) = filtfilt(b1, a1, freq_band_tmp);
            ca_pb_filt_new{ir}(:,i) = ca_pb_filt{ir}(round(len_filt/4)+1 : end-round(len_filt/4),i);% trimmed  ca
    end
end

% filter of 0~5Hz all
for ir = 1:trials
    freq_band_tmp = [fliplr(ca_pp_all{ir}(1:round(len_filt/4),:)); ca_pp_all{ir}; ...
        fliplr(ca_pp_all{ir}(end-round(len_filt/4)+1:end,:))];% extend signal for filtering
    ca_pb_filt_all{ir} = filtfilt(b1, a1, freq_band_tmp);
    ca_pb_filt_all_new{ir} = ca_pb_filt_all{ir}(round(len_filt/4)+1 : end-round(len_filt/4),:);% trimmed  ca
end
clear freq_band_temp
% filter of <0.05Hz power
for ir = 1:trials
    freq_band_tmp = [fliplr(ca_pp_bls_ave{ir}(:,1:round(len_filt/4)))'; ca_pp_bls_ave{ir}'; ...
        fliplr(ca_pp_bls_ave{ir}(:,end-round(len_filt/4)+1:end))'];% extend signal for filtering
    ca_pp_bls_filt{ir} = filtfilt(b1, a1, freq_band_tmp);
    ca_pp_bls_filt_new{ir} = ca_pp_bls_filt{ir}(round(len_filt/4)+1 : end-round(len_filt/4),:);% trimmed  ca
end
close all


%% Ca 1D plot

% ca time course(0.01~0.1Hz) 
for ir = 1:trials
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
for ir = 1:trials
h = figure;
    for i = 1:5
        subplot(5,1,i);
        plot(time2,zscore(ca_pp_ave{ir}(:,i)));
        xlim([0,640]); 

        title(['Ca power profile ',num2str(i-1),'~',num2str(i), ' (Hz)']);
        mean_power(:,i) = mean(ca_pp_ave{ir}(:,i),1);
        txt = ['\mu(dB) = ', num2str(mean_power(:,i))];
        legend(txt,'FontSize',12,'TextColor','red','FontSize',10);
        legend('boxoff');
    end
    saveas(h,['ca power profile trail ',num2str(scans(ir)),'.jpg']);
end
% Ca Power profile (0.01~0.1Hz)[separate]
for ir = 1:trials
    h = figure;
    for i = 1:5
        subplot(5,1,i);
        plot(time2,zscore(ca_pb_filt_new{ir}(:,i)));
        xlim([0,640]);
        ylabel('Power(normalized)');
        title(['Filtered Ca power profile ',num2str(i-1),'~',num2str(i), ' (Hz)']);
    end
    saveas(h,['Filtered Ca power profile ',num2str(scans(ir)),'.jpg']);
end
% Ca power profile (0.01~0.1Hz)[all averaged]
for ir = 1:trials
    h = figure;
    plot(time2,zscore(ca_pb_filt_all_new{ir}));
    xlim([0,640]);
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
for ir = 1:trials
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
for ir = 1:trials
    for is = 1:nslice
        freq_band_tmp = [fliplr(fmri_trim{ir}{is}(:,1:round(len_fmri/4)))'; fmri_trim{ir}{is}'; ...
        fliplr(fmri_trim{ir}{is}(:,end-round(len_fmri/4)+1:end))'];% extend signal for filtering
        fmri_vf{ir}{is} = filtfilt(b2, a2, freq_band_tmp);
        fmri_vf{ir}{is} = fmri_vf{ir}{is}(round(len_fmri/4)+1 : end-round(len_fmri/4),:)';
    end
end
% epoch-wise filtering
for ir = 1:trials
    for is = 1:nslice
        freq_band_tmp = [fliplr(per_fmri{ir}{is}(:,1:round(len_fmri/4)))'; per_fmri{ir}{is}'; ...
        fliplr(per_fmri{ir}{is}(:,end-round(len_fmri/4)+1:end))'];% extend signal for filtering
        fmri_pf{ir}{is} = filtfilt(b2, a2, freq_band_tmp);
        fmri_pf{ir}{is} = fmri_pf{ir}{is}(round(len_fmri/4)+1 : end-round(len_fmri/4),:)';
        for i = 1:nepoch
            ep_filt{ir}{is}(:,:,i) = fmri_pf{ir}{is}(:,1+(i-1)*nY:i*nY);
        end
        filt_ep{ir}{is} = mean(ep_filt{ir}{is},3);
    end
end
for ir = 1:trials
    h = figure;
    for is = 1:nslice
        subplot(3,1,is)
        
        for iv = nvoxel:-1:1
            plot(t_ep,100*filt_ep{ir}{is}(iv,:),'Color',clr_ep(iv,:));
            hold on;
            xlabel('times(s)');
            ytickformat('%g%%');
            ylabel('BOLD change');
        end
        h.Position = [1630,1330,400,700];
        title(['Epoch BOLD filt trial#',num2str(scans(ir)),' slice#',num2str(is)]);
    end
    saveas(h,['Epoch BOLD response trial#',num2str(scans(ir)),'.jpg']);
end

for ir = 1:trials
    for is = 1:nslice
        h = figure;
        imagesc(t_ep,[],100*filt_ep{ir}{is});
        colormap jet;
        colorbar;
        caxis([-1,3]);
        title(['Epoch BOLD 2D trial#',num2str(scans(ir)),' S#',num2str(is)]);
        saveas(h,['Epoch BOLD 2D trial#',num2str(scans(ir)),' S#',num2str(is),'.jpg']);
    end
end


% filtered fmri(0.01~0.1Hz)
for ir = 1:trials
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
%% xcorr(fMRI & Ca time-course)

    lg = 10;
    lg2 = 100;


for is = 1:nslice
      figure;
    for ir = 1:trials
        ca_down(ir,:) = downsample(ca_filt{ir},fs_ca*TR);
        [xcf(ir,is,:),lag_1] = xcorr(zscore(ca_down(ir,:)),zscore(fmri_filt{ir}{is}),lg2,'coeff');% calcium leading fmri by how long
        plot(lag_1*TR,squeeze(xcf(ir,is,:)));
        hold on;
        xlabel('Lag time(sec)');
        ylabel('correlation coefficient');
        title(['cross-correlation(fmri&ca) slice',num2str(is)]);
        [pks(ir,is,:),locs(ir,is,:)] = max(xcf(ir,is,:));
        locs(ir,is,:) = (locs(ir,is,:)-(lg2+1))*TR;
        plot(locs(ir,is,:),pks(ir,is,:),'*');
        xlim([-10,10])
    end
end


%% xcorr fmri slice
for ir = 1:trials
     h = figure;
     i = 0;
    for is = 1:nslice
        for si = 1:nslice
            i = i+1;
            sp_h = subplot(3,3,i);
            sp_h.Position = sp_h.Position+[0,0,0.03,0.03];
            
            [rho(ir,is,:,:),pval(ir,is,:,:)] = corr(fmri_vf{ir}{is}',fmri_vf{ir}{si}');
            imagesc(squeeze(rho(ir,is,:,:)));
            colormap jet;
            colorbar;
            caxis([0,1]);
            title(['Trials# ',num2str(scans(ir)),'slice#',num2str(is) ,'and',num2str(si) ]);
            
        end
    end
    h.Position = [300,300,1000,500];
    saveas(h,['Correlation of fMRI trials# ',num2str(scans(ir)),'.jpg']);
end
% cross slice fmri correlations
for ir = 1:trials
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
err_min = min(xcor_fmri(:,:,i),[],1);
err_max = max(xcor_fmri(:,:,i),[],1);
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

%% average Ca^2^+ frequency power-band
for ir = 1:trials
    ca_pb_filt_xc{ir} = mean(ca_pb_filt_new{ir},2);%1~5Hz [sperated]
end
%% xcorr (fMRI & Ca power-profile)
% xcorr (fmri & ca_power profile) [averaged] 
for is = 1:nslice
      h=figure;
      H = gca;
    for ir = 1:trials
        [xcf_test(ir,is,:),lag_2] = xcorr(zscore(ca_pb_filt_xc{ir}),zscore(squeeze(fmri_filt{ir}{is})),lg2,'coeff');% for 0-5Hz band
        lag = [-lg:lg/lg2:lg];
        h1= plot(lag,squeeze(xcf_test(ir,is,:)),'Color',[.65,.65,.65],'LineWidth',1.5);
        hold on;
        xlabel('Lag time(sec)');
        ylabel('correlation coefficient');
        title(['Averaged Ca^2^+(1~5Hz) & fMRI corr slice',num2str(is)]);
        [~,locs_pp(ir,is,:)] = max(abs(xcf_test(ir,is,:)));
        pks_pp(ir,is,:) = xcf_test(ir,is,locs_pp(ir,is,:));
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
%         legend;
        set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
saveas(h,['xcorr slice',num2str(is),'.jpg']);
end
save(['locs_pp'],'locs_pp');
save(['pks_pp'],'pks_pp');
% voxel-specific xcorr
% for ir = 1:trials    
%     for is = 1:nslice
%      test333{ir}{is} = downsample(fmri_vf{ir}{is}',1/TR)';
%     end
% end
%% Ca & fMRI & Ca power profile
for is = 1:nslice
    for ir = 1:trials
        h = figure;
        
        subplot(311);
        plot(time2,zscore(ca_pb_filt_all_new{ir}),'LineWidth',2,'Color',[0 0.4470 0.7410]);
        ylabel('Ca^2^+(1~5Hz)','fontweight','bold'); 
        xlim([0,640]);
        box off;
        title(['Lag time(s):',num2str(locs_pp(ir,is,:)),' CC:',num2str(pks_pp(ir,is,:))]);
        
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
        title(['Lag time(s):',num2str(locs(ir,is,:)),' CC:',num2str(pks(ir,is,:))]);
        xlim([0,640]);
        box off;
%         sgtitle(['Ultra-slow fluctuation in trial#',nun2str(scans(ir))]);

        saveas(h,['fmri_ca(1~5Hz) trail',num2str(scans(ir)),'slice',num2str(is),'.jpg']);
    end
    
end
close all


%% baseline , activity & BOLD comparison

cache_locs = [];
cache_locs_pp = [];
cache_pks = [];
cache_pks_pp = [];



for is = 1:nslice
    h=figure;
for ir = 1:trials
    
    
    if any(g_scans==scans(ir))
        plot(locs(ir,is,:),pks(ir,is,:),'k.','MarkerSize',10);
        hold on;
        plot(locs_pp(ir,is,:),pks_pp(ir,is,:),'r.','MarkerSize',10);
        line([locs(ir,is,:),locs_pp(ir,is,:)],[pks(ir,is,:),pks_pp(ir,is,:)],'Color','k');
        
        cache_locs= [cache_locs,locs(ir,is,:)];
        cache_pks = [cache_pks,pks(ir,is,:)];
        cache_locs_pp = [cache_locs_pp,locs_pp(ir,is,:)];
        cache_pks_pp = [cache_pks_pp,pks_pp(ir,is,:)];
    end
    legend('Baseline correlation','Power correlation');
    
end
title(['coupled comparison slice#',num2str(is)]);
xlim([-10,10]);
ylim([-1,1]);
saveas(h,['coupled comparison slice#',num2str(is),'.jpg']);
end

for is = 1:nslice
    h=figure;
for ir = 1:trials
    if ~any(g_scans==scans(ir))
        plot(locs(ir,is,:),pks(ir,is,:),'k.','MarkerSize',10);
        hold on;
        plot(locs_pp(ir,is,:),pks_pp(ir,is,:),'r.','MarkerSize',10);
        line([locs(ir,is,:),locs_pp(ir,is,:)],[pks(ir,is,:),pks_pp(ir,is,:)],'Color','k');
    end
end
title(['uncoupled comparison slice#',num2str(is)]);
xlim([-10,10]);
ylim([-1,1]);
saveas(h,['uncoupled comparison slice#',num2str(is),'.jpg']);
end
%%
for ir = 1:trials
    for is = 1:nslice
        for iv = 1:nvoxel
         [coeff_vo{ir}{is}(:,iv),lag_2] = xcorr(zscore(ca_pb_filt_xc{ir}),zscore(fmri_vf{ir}{is}(iv,:)'),lg2,'coeff');
        end
    end
end

% mean layer-specific lag time function handle
mean_lg = @(a,b) mean(a)-b; % the most frequent presence to be assigned
for ir = 1:trials
    for is = 1:nslice
            [~,lag_max_temp{ir}{is}] = max(abs(coeff_vo{ir}{is}),[],1);% GENERAL
            for iv = 1:nvoxel
            test111(iv) = coeff_vo{ir}{is}(lag_max_temp{ir}{is}(iv),iv);
            coeff_max(ir,is,:) = test111(iv);
            end
        for il = 1:5
            lag_max{ir}{is}(:,il) = round(mean_lg(lag_max_temp{ir}{is}(:,voxel_idx==il & abs(test111)>0.2),(lg2))/lg,1);
        end
        lag_max_tmp(ir,is,:) = lag_max{ir}{is};% pass cell to matrix
        lag_max_array(ir,is,:) = (lag_max_temp{ir}{is}-(lg2+1))/lg;
    end
end
% 1D xcorr [layer-specific]
Colour = colormap(hsv(5));
for ir = 1:trials
    for is = 1:nslice
        h = figure;
        xlim([-10,10]);


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
save(['l_coeff'],'coeff_vo');
% scatter layer-specific
for is = 1:nslice
    figure;
    for ir = 1:trials
        plot(1:5,lag_max{ir}{is}/lg,'*');
        xticks([1,2,3,4,5]);
        hold on
    end
end
for is = 1:nslice
    for i = 1:5
        lg_ave_t(is,i) = nanmean(lag_max_tmp(:,is,i),1);
        lg_std_t(is,i) = nanstd(lag_max_tmp(:,is,i),1);
    end
end
lg_ave = squeeze(lg_ave_t)';
lg_std = squeeze(lg_std_t)';
close all
%% error bar plot


%ALTERNATES
h = figure;
% x_group = [1:1:5];
% x_group = categorical({'L1','L2/3','L4','L5','L6'});
b = bar(lg_ave,.8,'FaceColor','flat','EdgeColor','flat');
hold on;

% Note: 'XOffset' was not documented
xCnt = (get(b(1),'XData') + cell2mat(get(b,'XOffset'))).';
er = errorbar(xCnt, lg_ave, lg_std/trials, 'k', 'LineStyle','none','LineWidth',3);

Fclr = [51, 153, 255;255, 51, 51;255, 153, 51]/255;
for ir = 1:trials
    for k = 1:3
        p(k) = plot(xCnt(:,k),squeeze(lag_max_tmp(ir,k,:)),'.','MarkerFaceColor',Fclr(k,:),'MarkerEdgeColor',Fclr(k,:));
        p(k).MarkerSize = 8;
%         p(k).MarkerFaceColor = Fclr(k,:);
        hold on;
    end
end
ylim([-10,10]);
legend({'Slice#1','Slice#2','Slice#3'});
legend('boxoff');
box off;
ax = gca;
ax.XTickLabel = {'L1','L2/3','L4','L5','L6'};
ax.YLabel.String = 'Lag time';
title('SEM bar of animal');
saveas(h,'SEM bar of animal.jpg');
%% 2D xcorr
for ir = 1:trials
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

for ir = 1:trials
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
if ave == 0
for is = 1:nslice
    h = figure;
    xcf_max = max(squeeze(xcf_test(:,is,:)),[],1);
    xcf_min = min(squeeze(xcf_test(:,is,:)),[],1);
    fill([lag';flipud(lag')],[xcf_min';flipud(xcf_max')],[0.3010 0.7450 0.9330],'LineStyle','none');
    line(lag,squeeze(mean(xcf_test(:,is,:),1)),'Color',[0,0,1],'LineWidth',2);
    [~,idx_1(is,:)] = max(mean(xcf_test(:,is,:),1));
    xline((idx_1(is,:)-lag_2(end))/lg,'LineWidth',1.5,'LineStyle','--','Color',[0 0.4470 0.7410]);
    xlabel('Lag time(s)');
    ylim([-1,1]);
    ylabel('Correlation coefficient');
    title(['Correlation of Ca^2^+(1~5Hz) power profile & fMRI slice# ',num2str(is)]);
    box off;
    legend('All trails','mean of all trails',['Averaged Peak Lag = ',num2str((idx_1(is,:)-lag_2(end))/lg),' sec']);
    legend('boxoff');
%     saveas(h,['Correlation of Ca(0.01~10Hz) power profile to fMRI slice# ',num2str(is),'.jpg']);
end
end

% fre-dependent xcorrelation (Ca & fMRI)

for i = 1:nfreq
  h = figure;
    for ir = 1:trials

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
        legend('Maximum lag');
        legend('boxoff');
%         legend;
    end
    text(0.8,9,'slice 1');
    title(['Ca^2^+(',num2str(i-1),'~',num2str(i),'Hz) power band & fMRI correlation']);
    saveas(h,['Ca^2^+(',num2str(i-1),'~',num2str(i),'Hz) power band & fMRI correlation.jpg']);
end
close all

clear  locs_pp_new

% xcorr (2~3Hz & <0.05Hz)
h = figure;
for ir =1 : trials
% for ir = 2
% figure;
    [coeff_1(:,ir),~] = xcorr(zscore(ca_pp_bls_filt_new{ir}),zscore(ca_pb_filt_all_new{ir}),lg2,'coeff');
    [pks_pb_new(ir,:),locs_pb(ir,:)] = max(coeff_1(:,ir));
    locs_pb(ir,:) = (locs_pb(ir,:)-(lg2+1))/lg;
    plot(lag,coeff_1(:,ir),'Color',[0.75 0.75 0.75]);
    hold on;
    plot(locs_pb(ir,:),pks_pb_new(ir,:),'r*');
    xlabel('lag time(s)');
    yline(0.2);
    ylabel('coeff');
    box off;
    title('xcorr between baseline & 1~5Hz Ca^2^+');
    ylim([-0.5,1]);
%     legend;
    title(num2str(scans(ir)));
end
text(2,0.4,['N = ',num2str(ir)],'FontSize',14);
saveas(h,'xcorr between baseline & 1~5Hz Ca^2^+ .jpg');
%% file I/O
temp_test = table(pks,locs,pks_pp,locs_pp);
temp_stat = table(lg_ave,lg_std);
temp_lr = table(coeff_max,lag_max_array);
temp_fmri = table(per_fmri,fmri_ep,fmri_max,onset);
temp_ca = table(per_ca,ca_ep,ca_match,ca_max);
save('all_lg_evoke.mat','temp_test'); 
save('stat_lg_evoke.mat','temp_stat');
save('layer_lag_evoke.mat','temp_lr');
save(['Ave_corr_evoke'],'xcf_test');
save('locs.mat','cache_locs');
save('locs_pp.mat','cache_locs_pp');
save('pks.mat','cache_pks');
save('pks_pp.mat','cache_pks_pp');
save('temp_fmri.mat','temp_fmri');
save('temp_ca.mat','temp_ca');
