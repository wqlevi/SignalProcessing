%% Initializing parameters
addpath('C:\Users\wangqi\Documents\Lab\Demo\Fid2plot') % the path where Brucker functions were stored
path = 'D:\27072020_yi\'; % targeted path of data 
addpath('C:\Users\wangqi\Documents\Lab\Demo\Calcium');% some functions imported
% path = 'D:\05232019\05232019_fmri_analyzed\20\fid';

params.complex = 1;
params.filter=[0,0,0,0];
params.phase = 0; 
params.zf = [0,0,0,0];
params.ft = [1,1,0,0];

cortex_depth = 2;% 2mm
spatial_res = 0.1; % see `method` 'Spatial_resolution'
len_cortex = cortex_depth/spatial_res; % cortex depth in samples
depth_cc = 5;% 0.5 mm of corpus callosum
cc_on = 0;% 0,w/o cc | 1,w/ cc
layer_idx = [1,2,2,2,2,3,3,4,4,4,4,4,4,5,5,5,5,5,5,5];
layer_name = {'L1','L2/3','L4','L5','L6'};

% sampling rates;
TR = 0.1;
Fs = 1/TR;
nepoch = 32;
nvoxels = 20;
% scans = [44,46,50]; % 14072020
scans = [34,35,36,38,39,40,41,42,43,44,45,46,47,48,49];%27072020
% scans = [31,32,33,34,35];%09102020
% scans = [40];
% scans = [48,49,56,57,58,59];
% scans = [63,64,65,66,67,68];
% scans = [32,34,35,36,37,38,39];
% scans = [34,35,36,38,39,40,41,42,43,44,45,46,47,48,49];
% scans = [61];
% scans = [54,55,56,57,58,59,60]; 
% scans = [35,40,81];% 16072020
% scans = [68];
% scans = [49,61,62,63]; %30072020_rs
% scans = [54,55,56,57,58,59,60]; % 30072020


prestim = 1/TR;

% scans = [77,78,80,81,82,83,84,90,91];% 16072020

ave = 1; % 1 | 0, 'ave' | 'rs'
rs = 0; % 1 | 0, 'rs' | 'evoked'
% I/O path:

mkdir([path,'\figs']);

save_path=[path,'\figs'];

%% Read raw fid
for ir = 1:length(scans)
        path_fmri = [path,num2str(scans(ir)),'\fid'];
        [res, p, pathstr] = BrReadImage_AutoRun(path_fmri, params);
        res = squeeze(res);
        fmri(ir,:,:,:,:) = squeeze(res); % trials*128*32*nslice*200
end
%% Average trails(evoked only)
if ave ==1
fmri = squeeze(mean(fmri,1));
disp(['Number of averaged trials: ',num2str(ir)])
end

%% Formating fid
if ave == 1
[nX, nY, ndur] = size(squeeze(fmri(:,:,1,:)));   
nslice = size(fmri,3);
trials = 1;
else
[nX, nY, ndur] = size(squeeze(fmri(1,:,:,1,:))); % nX = cortical depth | nY = epochs | nSlice = duration of epoch
nslice = size(fmri,4);
trials = length(scans);
end

% nslice = size(fmri,3);% ave ,size(,3)
for ir = 1:trials
    for is = 1:nslice
        if ave==1
        fmri_slice{ir}{is} = squeeze(fmri(:,:,is,:)); %128*32*200
        else  
        fmri_slice{ir}{is} = squeeze(fmri(ir,:,:,is,:));% ave
        end
        slice_fftc{is} = fft2c(fmri_slice{ir}{is});
        slice_ifft{is} = ifft2c(slice_fftc{is});
        fmri_slice{ir}{is} = fftc(slice_ifft{is},2);
       for k = 1:32
        for j = 1:200
            fmri_cm(:,k+32*(j-1)) = abs(fmri_slice{ir}{is}(:,k,j));
        end
       end
    fmri_tc{ir}{is} = fmri_cm - min(min(fmri_cm)); 
    
    end
end

t_fmri = 0:TR:(nY*ndur-1)*TR; % t_fmri initialization
% tSNR plot
for ir = 1:trials
    for is = 1:nslice
        h = figure; 
        set(h, 'WindowStyle', 'docked', 'NumberTitle', 'off', 'Name', ['tSNR '])
        plot(linspace(1,nX,nX),abs(mean(fmri_tc{ir}{is},2)));
        xlabel('cortical depth(mm)');
        xlim([1,128]);
        ylabel('intensity');
        title(['Cortical surface trail#',num2str(scans(ir))]);
    end
end
for ir = 1:trials
    for is = 1:nslice
        h = figure;
        plot(linspace(1,nX,nX),tSNR(abs(mean(fmri_tc{ir}{is},2))));
        ylabel('tSNR');
        xlabel('Cortical depth');
    end
end
%% Cutting off cortical depth
cd(save_path)
% CC border still untargeted
for ir = 1:trials
    for is = 1:nslice
        ave_tsnr(:,ir,is) = abs(mean(fmri_tc{ir}{is},2));
        ave_tsnr(:,ir,is) = ave_tsnr(:,ir,is) - min(ave_tsnr(:,ir,is));
%         half_idx(:,ir,is) = find(ave_tsnr(:,ir,is)>(max(ave_tsnr(:,ir,is))/2),1);
        half_idx(:,ir,is) = 39; % MANUALLY DEFINED CORTEX SURFACE
        half_idx_bot(:,ir,is) = half_idx(:,ir,is)+len_cortex; % bottom of cortex
        
        cc_idx_bot(:,ir,is) = half_idx_bot(:,ir,is)+depth_cc; % bottom of corpus callosum
    end
end

for ir = 1:trials
    for is = 1:nslice
        h = figure();
        subplot(1,2,1)
        imagesc(t_fmri,[],fmri_tc{ir}{is});
        colormap bone;
        hold on;
        yline(half_idx(:,ir,is),'g--','LineWidth',2);
        yline(half_idx_bot(:,ir,is),'r--','LineWidth',2);
        yline(cc_idx_bot(:,ir,is),'--','LineWidth',2);
        ylabel('cortical depth');
        xlabel('time(s)');
        title(['cortical map trail# ',num2str(ir),'slice#',num2str(is)]);
        subplot(1,2,2)
        plot(linspace(1,nX,nX),ave_tsnr(:,ir,is));
        hold on;
        l(:,1)=xline(half_idx(:,ir,is),'g--','LineWidth',3);
        l(:,2)=xline(half_idx_bot(:,ir,is),'r--','LineWidth',3);
        l(:,3)=xline(cc_idx_bot(:,ir,is),'--','LineWidth',3);
        xlabel('cortical depth');
        ylabel('time(s)');
        xlim([1,128]);
        set(gcf,'Position',[500 500 980 300 ]);
        title('cortical surface trimming');
        legend([l(:,1),l(:,2),l(:,3)],'Start cortex','end cortex','end CC')

        saveas(h,['Cutting off trail#',num2str(scans(ir)),'slice#',num2str(is),'.jpg']);
    end
end
% trim along cortical depth here
for ir = 1:trials
    for is = 1:nslice
        fmri_cm_cut{ir}{is} = fmri_tc{ir}{is}(half_idx(:,ir,is):half_idx_bot(:,ir,is)-1,:);
    end
    % trimmed depth including corpus callosum
    for is = 1:nslice
        fmri_cc_cut{ir}{is} = fmri_tc{ir}{is}(half_idx(:,ir,is):cc_idx_bot(:,ir,is)-1,:);
    end
end


%% tSNR

for is = 1:nslice
    tsnr_temp(:,ir,is) = squeeze(mean(fmri_cm_cut{ir}{is},2) ./ std(fmri_cm_cut{ir}{is},0,2));% computation: mean/stdev
end

for ir = 1:trials
    for is = 1:nslice
        figure(),
        plot(1:size(tsnr_temp,1),squeeze(tsnr_temp(:,ir,is)));
        xlim([1,size(tsnr_temp,1)]);
        title('tSNR of cortex');
    end
end
%% cut off [ABORTED]
for ir = 1:trials
    for is = 1:nslice
        norm_cm{ir}{is} = fmri_cm_cut{ir}{is} ./ max(max(fmri_cm_cut{ir}{is}));
        norm_cc{ir}{is} = fmri_cc_cut{ir}{is} ./ max(max(fmri_cc_cut{ir}{is})); % corpus callosum included
%         figure(),
%         cmin = min(min(norm_cm{ir}{is}));
%         cmax = max(max(norm_cm{ir}{is}));
%         imagesc(t_fmri,[],squeeze(norm_cm{ir}{is}));
%         xlabel('time(s)');
%         ylabel('cortical voxels');
%         colormap(jet);
%         caxis([cmin,cmax]);
%         set(gcf,'Position',[500,500,600,300]);
    end
end

%% Normalize fMRI baseline to 100%
disp('Normalizing baseline of timecourse to 100')
for ir = 1:trials
    for is = 1:nslice
        fmri_cm_de{ir}{is} = (fmri_cm_cut{ir}{is}./mean(fmri_cm_cut{ir}{is},2))*100;
    end
end
%% Epoch-wise BOLD
t_epoch = [0:TR:(ndur-1)*TR];
for ir = 1:trials
    for is = 1:nslice
        for i = 1:nepoch
            for iv  = 1:nvoxels
                prestim_ref{ir}{is}(iv,:,i) = mean(fmri_cm_de{ir}{is}(iv,1+(i-1)*ndur:prestim+(i-1)*ndur),2);
                fmri_cm_temp(iv,:,i) = (fmri_cm_de{ir}{is}(iv,1+(i-1)*200:i*200)-prestim_ref{ir}{is}(iv,:,i))./prestim_ref{ir}{is}(iv,:,i); % Epoch-wise: \Delta F/ F
%                 disp(prestim_ref)
            end
        end
        fmri_time_temp{ir}{is} = reshape(fmri_cm_temp,nvoxels,ndur*nepoch);
        fmri_cm_epoch{ir}{is} = mean(fmri_cm_temp,3); % epoch-ave, cortical indexed BOLD(Unit:100%)
    end
end

% plot the cortical epoch response
% 2D
for ir = 1:trials
for is = 1:nslice
    h = figure;
    imagesc(t_epoch,[],fmri_cm_epoch{ir}{is});
    colormap jet;
    colorbar;
    caxis([-(2e-3),8e-3]);
    saveas(h,['Epoch BOLD map trial#',num2str(scans(ir)), 'slice#',num2str(is),'.jpg']);
end
end
% 1D
COLOR = flip(colormap(jet(nvoxels)),1);
for ir = 1:trials
    for is = 1:nslice
        h = figure;
        for iv = 1:nvoxels
            v = plot(0:TR:(ndur-1)*TR,fmri_cm_epoch{ir}{is}(iv,:)*100,'Color',COLOR(iv,:));
            hold on;
        end
        xlabel('Time(s)');
        ylabel('Percentage(%)');
        ytickformat('percentage');
        title('Epoch-wise averaged BOLD response');
        colormap(COLOR);
        colorbar('Ticks',[0,1],'Direction','reverse','TickLabels',{'Surface','Bottom'});
        saveas(h,['Epoch-wise averaged BOLD response slice#',num2str(is),'.jpg'])
    end
end
%% Filtering(for 2D BOLD map)

order    = 4096; 
Fc1  = 0.01;   %low cut frequency in Hz 0.01 
Fc2 = 0.1;   %intend 0.1 Hz but consider transient band
N = order;
beta = 0.005;
win = kaiser(N+1, beta); %using kaiser window
flag = 'scale';  % Sampling Flag

% Calculate the coefficients using the FIR1 function.
b  = fir1(N, [Fc1, Fc2]/(Fs/2), 'bandpass', win, flag);

fvtool(b,1,'Fs',Fs)
Hd = dfilt.dffir(b);

% apply filters
switch cc_on
    case 0
        cortical_depth_map = norm_cm;
        save_name_2 = 'temperal-spatial fmri map';
        disp('filtering without Corpus callosum...')
    case 1
        cortical_depth_map = norm_cc;
        save_name_2 = 'temperal-spatial fmri map(with CC.)';
        disp('filtering with Corpus callosum...')
end 
if ave == 1
    save_name_3 = ' ave';
else 
    save_name_3 = ' ';
end
% dimension of vairables, dependent of signal depth
filtered_cortical_depth_map = zeros(size(cortical_depth_map{ir}{is}));
nVoxel = size(cortical_depth_map{ir}{is},1); % num of voxels
nTime = size(filtered_cortical_depth_map,2); % num of time points
cortical_depth_map_all = zeros(nVoxel,nTime,size(fmri,3));
for ir = 1:trials
    for is = 1:nslice
    %     cortical_depth_map = norm_cm{i};

        for voxel_num = 1 : nVoxel

            %prevent the signal fluctuation in front and behind part after filtering
            temp_filtered_cortical_depth_map1(voxel_num,:)=fliplr((cortical_depth_map{ir}{is}(voxel_num,2:N+1)));  % maintain continuity in level and slope
            temp_filtered_cortical_depth_map2(voxel_num,:)=fliplr((cortical_depth_map{ir}{is}(voxel_num,end-N:end-1)));

            filtered_cortical_depth_map_new(voxel_num,1:N)                 = temp_filtered_cortical_depth_map1(voxel_num,1:N);
            filtered_cortical_depth_map_new(voxel_num,N+1:nTime+N)         = (cortical_depth_map{ir}{is}(voxel_num,1:nTime));
            filtered_cortical_depth_map_new(voxel_num,nTime+N+1:nTime+2*N) = temp_filtered_cortical_depth_map2(voxel_num,1:N);

            filtered_cortical_depth_map_new(voxel_num,:) =  filter(Hd,squeeze(filtered_cortical_depth_map_new(voxel_num,:)));
            delay = mean(grpdelay(Hd));
            filtered_cortical_depth_map_new(voxel_num,:) = circshift(filtered_cortical_depth_map_new(voxel_num,:),(-1)*delay,2);
        end
        filtered_cortical_depth_map(:,1:nTime) =  filtered_cortical_depth_map_new(:,N+1:N+nTime);
        cortical_depth_map_all(:,:,ir,is) = abs(filtered_cortical_depth_map);
    end
end
clear temp_filtered_cortical_depth_map1 temp_filtered_cortical_depth_map2 filtered_cortical_depth_map_new filtered_cortical_depth_map

%% epoch-wise normalized colormap(cortical temporal BOLD map)
% cortical_depth_map = temp_cm_cut(:,:,3);


clear cortical_depth_map;
cortical_depth_map = cortical_depth_map_all;
nY =32;
pre_stim = 10;
for ir = 1:trials
    for is = 1:nslice
        for nRepeat = 1 : nY %32

            SI_base_mean = mean(abs(cortical_depth_map(:,1+ndur*(nRepeat-1):pre_stim+ndur*(nRepeat-1),ir,is)),2); % 1s off
            SI_i = abs(cortical_depth_map(:,1+ndur*(nRepeat-1):ndur*(nRepeat),ir,is));% single epoch(1/32) response colormap
            for t = 1 : size(SI_i,2)
                SI_percentage(:,t) = (SI_i(:,t) - SI_base_mean(:,1))./(SI_base_mean(:,1));% percentage change within single epoch 
            end
            Threshold_SI_MAP(:,1+ndur*(nRepeat-1):ndur*(nRepeat)) = SI_percentage;
        end

        norm_cortical_map(:,:,ir,is) = abs(Threshold_SI_MAP(:,:))./max(max(abs(Threshold_SI_MAP(:,:))));

        norm_cortical_map_temp(:,:,ir,is) = abs(cortical_depth_map(:,:,ir,is))./max(max(abs(cortical_depth_map(:,:,ir,is))));
        Threshold_SI_MAP = abs(cortical_depth_map(:,:,ir,is))./mean(norm_cortical_map_temp(:,:,ir,is),2);
        norm_cortical_map(:,:,ir,is) = abs(Threshold_SI_MAP(:,:))./max(max(abs(Threshold_SI_MAP(:,:))));
    end
end
% plot colormap of BOLD 
TA = 640;

for ir = 1:trials
    for is = 1:nslice
        h = figure; set(h, 'WindowStyle', 'docked', 'NumberTitle', 'off', 'Name', ['fMRI time course trail# ',num2str(ir),' slice ', num2str(is)])
        imagesc(t_fmri,[],norm_cortical_map(:,:,ir,is));
        cmax = max(max(norm_cortical_map(:,:,ir,is)));
        cmin = min(min(norm_cortical_map(:,:,ir,is)));
        colormap(jet);
        colorbar;
        xlim([0,t_fmri(end)]);
        xlabel('time(s)');
        ylabel('voxels');
        title([save_name_2,save_name_3]);
        caxis([0.2,cmax])
        set(gcf,'Position',[500 500 980 300 ]);
        saveas(h,[save_name_2,save_name_3,'trail# ',num2str(scans(ir)),' slice# ',num2str(is),'.jpg']);
    end
end
% clear SI_percentage Threshold_SI_MAP norm_cortical_map norm_cortical_map_temp
%% Time course of ave fmri
switch cc_on
    case 1 
        title_name = 'fmri time course including cc';
        title_name_filted = 'fmri time course including CC(0.01-1Hz)';
        fmri_timecourse = fmri_cc_cut;

    case 0
        title_name = 'fmri time course';
        title_name_filted = 'fmri time course(0.01-0.1Hz)';
        fmri_timecourse = fmri_cm_cut;

end
switch ave
    case 1
        title_name2 = ' averaged ';
    case 0
        title_name2 = ' single ';
end
        
for ir = 1:trials 
    if rs == 1
        break
    end
    for is = 1:nslice
    h = figure; 
    % set(h, 'WindowStyle', 'docked', 'NumberTitle', 'off', 'Name', ['fMRI time course slice', num2str(i)])
    plot(t_fmri,100*mean(fmri_time_temp{ir}{is},1),'k');
    xlim([0,t_fmri(end)]);
    ylabel('fMRI Percentage(%)');
    ytickformat('percentage');
    % xticks(0:6400/4*(1/TR):6400*(1/TR));
    xlabel('time(s)');
    % xticklabels({'0 sec','160 sec','320 sec','480 sec','640 sec'})

    title([title_name,title_name2,' #trail ',num2str(scans),' slice#',num2str(is)]);
    set(gcf,'Position',[500 500 980 300 ]);
    box off;
    saveas(h,[title_name,title_name2,' trail#',num2str(scans),' slice#',num2str(is),'.jpg']);
    end
end
for ir = 1:trials
    for is = 1:nslice
    mean_slice_1(ir,is,:) = squeeze(mean(abs(fmri_timecourse{ir}{is}),1));
    mean_slice_2(ir,is,:) = detrend(squeeze(mean_slice_1(ir,is,:)));
    end
end

% filtering fmri
bpf_low = 0.01;
bpf_high = 0.1;
[b,a] = butter(1,[bpf_low bpf_high]/(Fs/2),'bandpass');
len = size(fmri_timecourse{ir}{is},2);
for ir = 1:trials
    for is = 1:nslice
        fmri_ave{ir}{is} = mean(fmri_timecourse{ir}{is},1);
        temp = [fliplr(fmri_ave{ir}{is}(:,1:len/4)),fmri_ave{ir}{is},fliplr(fmri_ave{ir}{is}(:,end-len/4+1:end))];
        temp_filt1 = filtfilt(b,a,temp);
        fmri_ave_f{ir}{is} = temp_filt1(:,len/4+1:end-len/4);
    end
end
% epoch-wise filtering
len_ep = size(ep_norm1{ir}{is},2);
for ir = 1:trials
    for is = 1:nslice
        for iv  = 1:nvoxels
            fmri_ep_tmp{ir}{is}(iv,:) = fmri_time_temp{ir}{is}(iv,:);
            temp = [fliplr(fmri_ep_tmp{ir}{is}(iv,1:len/4)),fmri_ep_tmp{ir}{is}(iv,:),fliplr(fmri_ep_tmp{ir}{is}(iv,end-len/4+1:end))];
            temp_filt3 = filtfilt(b,a,temp);
            fmri_ep{ir}{is}(iv,:) = temp_filt3(:,len/4+1:end-len/4);
        end
    end
end

for ir = 1:trials
    for is = 1:nslice
        h = figure;
        for il = 1:5
            subplot(3,2,il)
            plot(t_epoch,100*mean(fmri_ep{ir}{is}(layer_idx == il,:),1));
            ylabel({layer_name{il},'Percentage(%)'});
        end
    end
end

for ir = 1:trials
    for is = 1:nslice
        h = figure;
        Clr = colormap (winter(3));
        plot(t_fmri,zscore(squeeze(mean_slice_2(ir,is,:))),'k');
        hold on;
        plot(t_fmri,zscore(fmri_ave_f{ir}{is}),'Color',Clr(3,:),'LineWidth',2);
        xlim([0,t_fmri(end)]);
        ylim([-4,4]);
        xlabel('time(s)');
        ylabel('fMRI a.u.');
        grid on;
        box off;
        title([title_name_filted,title_name2,'trail#',num2str(scans(ir)),' slice#',num2str(is)]);
        set(gcf,'Position',[500 500 980 300 ]);
        saveas(h,[title_name_filted,title_name2,'trail#',num2str(scans(ir)),' slice#',num2str(is),'.jpg']);
    end
end
% spectral analysis

   
%% PSD
% make error bar plot  
for is = 1:nslice

    for ir = 1:trials
        test1(ir,is,:) = squeeze(mean_slice_1(ir,is,:));
        test1(ir,is,:) = test1(ir,is,:) - mean(test1(ir,is,:));
        [pxx(ir,is,:),f] = pwelch(squeeze(test1(ir,is,:)),hamming(size(test1,3)),0,[],1/TR);
    end
end
bsl = (mean(pxx(ir,is,0.9<f<1)));
pxx_norm = pxx/bsl;
for ir = 1:trials
        h = figure;
        clr = {'b','r'};
    for is = 1:2
%         plot(f(f>0&f<0.12),10*log10(pxx(f>0&f<0.12,:)),'Color','k');
        plot(f(f<0.5),squeeze(pxx_norm(:,is,f<0.5)),'Color',clr{is});
        ax = gca();
        ax.XScale = 'log';
        xticks([0,0.01,0.02,0.05,0.1,0.16,0.2,0.3,0.4]);
        xlim([f(2),f(size(f(f<0.5),1))]);
        hold on;
%         h.Position = [500 500 980 300 ];
    end
        xlabel('Frequency(Hz)');
        ylabel('Norm power');
        title(['spectrum of trail#',num2str(ir),' slice#',num2str(is)]);
        saveas(h,['spectrum of ',num2str(scans(ir)),'slice#',num2str(is),'.jpg']);
end

%         title(['spectrum of all slice#',num2str(is)]);
   
% log scaled PSD
for ir = 1:trials
        h = figure;
        clr = {'b','r'};
    for is = 1:2
%         plot(f(f>0&f<0.12),10*log10(pxx(f>0&f<0.12,:)),'Color','k');
        plot(f(f<0.5),10*log10(squeeze(pxx_norm(:,is,f<0.5))),'Color',clr{is});
        ax = gca();
        ax.XScale = 'log';
        xticks([0,0.01,0.02,0.05,0.1,0.16,0.2,0.3,0.4]);
        xlim([f(2),f(size(f(f<0.5),1))]);
        hold on;
%         h.Position = [500 500 980 300 ];
    end
        xlabel('Frequency(Hz)');
        ylabel('Norm power');
        title(['spectrum of trail#',num2str(ir),' slice#',num2str(is)]);
end


% layer-specific PSD

for ir = 1:trials
    for is = 1:nslice
        fmri_layer{ir}{is}(1,:) = mean(fmri_timecourse{ir}{is}(1,:),1);
        fmri_layer{ir}{is}(2,:) = mean(fmri_timecourse{ir}{is}(2:5,:),1);
        fmri_layer{ir}{is}(3,:) = mean(fmri_timecourse{ir}{is}(6:7,:),1);
        fmri_layer{ir}{is}(4,:) = mean(fmri_timecourse{ir}{is}(8:13,:),1);
        fmri_layer{ir}{is}(5,:) = mean(fmri_timecourse{ir}{is}(14:19,:),1);
    end
end
for ir = 1:trials
    for is = 1:nslice
        h = figure;
        for i = 1:5
            subplot(3,2,i);
            plot(t_epoch,100*mean(fmri_cm_epoch{ir}{is}(layer_idx == i,:),1),'Color',clr(i,:));
            ylabel({'Percentage (%)',layer_name{i}});
        %     ylim([]);
             xlabel('time(s)');
        end
        saveas(h,['Epeoch timecourse slice#',num2str(is),'.jpg']);
    end
end


clr = colormap(hsv(5));
for is = 1:nslice
   h = figure;
for ir = 1:trials
    for i = 1:5
        subplot(2,3,i)
        test = fmri_layer{ir}{is}(i,:);
        test = test - mean(test);
        [pxx,f] = pwelch(test,hamming(size(test,2)/4),0,[],1/TR);
        plot(f(f>0&f<0.1),pxx(f>0&f<0.1,:),'Color',clr(i,:));
        hold on;
        xlim([f(2),0.12]);
%         ylim([0,10])
        text(0.08,4,['N = ',num2str(length(scans))]);
        xlabel('frequency(Hz)');
        ylabel('intensity (dB/Hz)');
        title(['spectrum of trail#',num2str(ir),' slice#',num2str(is)]);
%         title(['spectrum of all slice#',num2str(is)]);
    end
end
%     saveas(h,['spectrum of slice#',num2str(is),'.jpg']);
end
%% Epoch-wise BOLD[raw fmri]
for ir = 1:trials
    for is = 1:nslice
        for i = 1:nepoch
            ep_temp1(:,:,i) = fmri_timecourse{ir}{is}(:,1+(i-1)*200:i*200);
            for iv = 1:nvoxels
                ep_temp1(iv,:,i) =( ep_temp1(iv,:,i) - mean(ep_temp1(iv,1:prestim,i),2))./mean(ep_temp1(iv,1:prestim,i),2);
            end
        end
        ep_norm1{ir}{is} = mean(ep_temp1,3);
    end
end
% for ir = 1:trials
%     for is = 1:nslice
%         h = figure;
%         for iv = 1:nvoxels
%         plot(0:TR:(ndur-1)*TR,ep_norm1{ir}{is}(iv,:),'Color',COLOR(iv,:));
%         xlabel('time(s)');
%         ylabel('\Delta F/F');
%         title(['raw epoch-wise dynamic trail#',num2str(scans(ir)),'slice#',num2str(is),'.jpg']);
%         hold on;
%         end
% %         saveas(h,['raw epoch-wise dynamic trail#',num2str(scans(ir)),'slice#',num2str(is),'.jpg']);
%     end
% end


for ir = 1:trials
   for is = 1:nslice
    h = figure;
    imagesc(ep_norm1{ir}{is});
    xlabel('time(s)');
    ylabel('cortical depth(voxels)');
    title(['Epoch-wise raw response trail# ',num2str(scans(ir)),' slice# ',num2str(is)]);
   saveas(h,['raw reponse epoch slice# ',num2str(is),'.jpg']);
   end
end
%% Epoch-wise BOLD[filtered fmri]
for ir = 1:trials
    for is = 1:nslice
        for i = 1:nepoch
            ep_temp(:,:,i) = norm_cortical_map(:,1+(i-1)*200:i*200,ir,is);
            for iv = 1:nvoxels
                ep_temp(iv,:,i) = (ep_temp(iv,:,i) - mean(ep_temp(iv,1:prestim,i),2))./mean(ep_temp(iv,1:prestim,i),2);
            end
        end
        ep_norm{ir}{is} = mean(ep_temp,3);
    end
end
for ir = 1:trials
    for is = 1:nslice
        h = figure;
        for iv = 1:nvoxels
            sbp_h = subplot(5,4,iv);
            plot(0:TR:(ndur-1)*TR,ep_norm{ir}{is}(iv,:),'k');
%             sbp_h.Position = sbp_h.Position + [0,0,0,0.03];
            title(['Voxel# ',num2str(iv)]);
            if mod(iv-1,4)
                ylabel(' ');
            else
                ylabel('\Delta F/F');
            end

            if iv<=nvoxels-4
                xlabel(' ');                
            else
                xlabel('Time(s)');
            end

        end
        xlabel('time(s)');
%         ylabel('\Delta F/F');
%         sgtitle(['filtered epoch-wise dynamic trail#',num2str(scans(ir)),'slice#',num2str(is)]);
        h.Position = [0,0,900,650];
%         legend;
        saveas(h,['filtered epoch-wise dynamic trail#',num2str(scans(ir)),'slice#',num2str(is),'.jpg']);
    end
end

for ir = 1:trials
   for is = 1:nslice
    h = figure;
    imagesc(ep_norm{ir}{is});
    colormap(jet);
    caxis([-.01,.05])
    xlabel('time(s)');
    ylabel('cortical depth(voxels)');
    colorbar;
    title(['Epoch-wise filtered response trail# ',num2str(scans(ir)),' slice# ',num2str(is)]);
   saveas(h,['Filtered reponse epoch #slice ',num2str(is),'.jpg']);
   end
end

%% Error bar plot[SNR]
% Manually type in intensity of each ROI of all animals:
% % 25 min scans
% roi(:,1) = [164,191];
% roi(:,2) = [24.5,71.4];
% roi(:,3) = [4.7,2.84];
% 6 min scans
roi(:,1) = [25.5,49.1,37.2]; % roi 2 in anatomic
roi(:,2) = [175,165,111]; % roi 1 in anatomic
roi(:,3) = [5.12,2.7,2.81]; % background signals(as reference) 
% function handle of SNR
snr = @(s,b) (s)/sqrt(b);

% calculation of SNR
for j = 1:2% 2 rois
    for i = 1:size(roi,1) % trails 
        snr_re(i,j) = snr(roi(i,j),roi(i,3));
    end
end
for j = 1:2
    ave(:,j) = mean(snr_re(:,j));
    snr_std(:,j) = std(snr_re(:,j));
end


% componds of plots
h = figure;
x_group = [1:1:2];% two groups with 1 interval
b = bar(x_group,ave,.3,'FaceColor','flat','EdgeColor','flat');
b.CData(1,:) = [0.8500 0.3250 0.0980];%clr bar#1
b.CData(2,:) = [0 0.4470 0.7410];%clr bar#2

hold on;
er = errorbar(x_group,ave,snr_std);
er.Color = [0,0,0];
er.LineStyle = 'none';
er.LineWidth = 2;
for k = 1:3
%     p = plot(1:2,snr_re(k,:),'*');
    p = plot(1,snr_re(k,1),'.',2,snr_re(k,2),'.');
    p(1).Color = [206,132,103]/255;
    p(2).Color = [131,177,210]/255;
    p(1).MarkerSize = 25;
    p(2).MarkerSize = 25;
    hold on;
end
hold off;
ax = gca;
ax.XTickLabel = {'ROI 1','ROI 2'};
ax.YLabel.String = 'SNR';
box off;
saveas(h,'snr plot.jpg');



