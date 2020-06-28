%% Initializing parameters
addpath('C:\Users\wangqi\Documents\Lab\Demo\Fid2plot') % the path where Brucker functions were stored
path = 'C:\Users\wangqi\Downloads\70\fid'; % targeted path of data 
% path = 'D:\05232019\05232019_fmri_analyzed\20\fid';

params.complex = 1;
params.filter=[0,0,0,0];
params.phase = 0; 
params.zf = [0,0,0,0];
params.ft = [1,1,0,0];

cortex_depth = 2;% 2mm
spatial_res = 0.1; % see `method` 'Spatial_resolution'
len_cortex = 2/spatial_res; % cortex depth in samples
depth_cc = 5;% 0.5 mm of corpus callosum

% sampling rates;
TR = 0.1;
Fs = 1/TR;

%% Read raw fid
[res, p, pathstr] = BrReadImage_AutoRun(path, params);
%% Formating fid
fmri = squeeze(res);
for i = 1:size(fmri,3)% slices
%     slice{i} = abs(fmri(:,:,i,:));
    fmri_slice{i} = squeeze(fmri(:,:,i,:));
    slice_fftc{i} = fft2c(fmri_slice{i});
    slice_ifft{i} = ifft2c(slice_fftc{i});
    fmri_slice{i} = fftc(slice_ifft{i},2);
%     [why not?] how does the `reshape` function organize matrix?
%     fmri_slice{i} = reshape(slice_fftc{i},size(slice_fftc{i},1),(size(slice_fftc{i},2)*size(slice_fftc{i},3))); % [why?]by doing fft, the profile is more salient
    for k = 1:32
        for j = 1:200
            fmri_cm(:,k+32*(j-1)) = abs(fmri_slice{i}(:,k,j));
        end
    end
fmri_tc{i} = fmri_cm; 
[nX, nY, nSlice] = size(squeeze(fmri(:,:,i,:)));
end
% tSNR plot
for i = 1:size(fmri,3)
    h = figure; set(h, 'WindowStyle', 'docked', 'NumberTitle', 'off', 'Name', ['tSNR '])
    plot(linspace(1,size(fmri_tc{1},1),size(fmri_tc{1},1)),abs(mean(fmri_tc{i},2)));
    xlabel('cortical depth(mm)');
    xlim([0,128]);
    ylabel('intensity');
    title('Cortical surface');
end

%% cortical mapping [aborted]
% for i = 1:size(fmri,3)
% [nX, nY, nSlice] = size(squeeze(fmri(:,:,i,:)));
% % temp_img = squeeze(fmri(:,:,i,:));
% % selected_ksp = fft2c(temp_img);
% % temp_lp = ifft2c(selected_ksp);
% % t0_cm = fftc(temp_lp,2);
% for k = 1:32
%     for j = 1:200
%         temp_cm(:,k+32*(j-1)) = abs(t0_cm(:,k,j));
%     end
% end
% temp_cm_new(:,:,i) = temp_cm; 
% end
%% Cutting off cortical depth
% CC border still untargeted

for i = 1:size(fmri,3)
    ave_tsnr(:,i) = abs(mean(fmri_tc{i},2));
    ave_tsnr(:,i) = ave_tsnr(:,i) - min(ave_tsnr(:,i));
    half_idx(:,i) = find(ave_tsnr(:,i)>(max(ave_tsnr(:,i))/2),1);
    half_idx_bot(:,i) = half_idx(:,i)+len_cortex;
    cc_idx_bot(:,i) = half_idx_bot(:,i)+depth_cc;
end

for i = 1:size(fmri,3)
    figure(),
    subplot(1,2,1)
    imagesc(fmri_tc{i});
%     imagesc(abs(fmri_slice{i}));
    colormap bone;
    hold on;
    yline(half_idx(:,i),'g--','LineWidth',2);
    yline(half_idx_bot(:,i),'r--','LineWidth',2);
    yline(cc_idx_bot(:,i),'--','LineWidth',2);
    ylabel('cortical depth');
    xlabel('time(s)');
    subplot(1,2,2)
    plot(linspace(1,size(ave_tsnr(:,i),1),size(ave_tsnr(:,i),1)),ave_tsnr(:,i));
    hold on;
    l(:,1)=xline(half_idx(:,i),'g--','LineWidth',3);
    l(:,2)=xline(half_idx_bot(:,i),'r--','LineWidth',3);
    l(:,3)=xline(cc_idx_bot(:,i),'--','LineWidth',3);
    xlabel('cortical depth');
    ylabel('time(s)');
    xlim([0,128]);
    set(gcf,'Position',[500 500 980 300 ]);
    title('cortical surface trimming');
    legend([l(:,1),l(:,2),l(:,3)],'Start cortex','end cortex','end CC')
end
% trim along cortical depth here
for i = 1:size(fmri,3)
    temp_cm_cut{i} = fmri_tc{i}(half_idx:half_idx_bot-1,:);
end
% trimmed img including corpus callosum
for i = 1:size(fmri,3)
    temp_cc_cut{i} = fmri_tc{i}(half_idx:cc_idx_bot,:);
end


%% tSNR

for i = 1:size(fmri,3)
    tsnr_temp(:,i) = squeeze(mean(temp_cm_cut{i},2) ./ std(temp_cm_cut{i},0,2));% computation: mean/stdev
end

for i = 1:size(fmri,3)
    figure(),
    plot(1:size(temp_cm_cut{i},1),tsnr_temp(:,i));
    xlim([1,size(temp_cm_cut{i},1)]);
    title('tSNR');
end
%% cut off 
for i = 1:size(fmri,3)
    norm_cm{i} = temp_cm_cut{i} ./ max(max(temp_cm_cut{i}));
    norm_cc{i} = temp_cc_cut{i} ./ max(max(temp_cc_cut{i})); % corpus callosum included
    figure(),
    cmin = min(min(norm_cm{i}));
    cmax = max(max(norm_cm{i}));
    imagesc(squeeze(norm_cm{i}));
    colormap(jet);
    caxis([cmin,cmax]);
end
%% Filtering

order    = 4096; %use it now

Fc1  = 0.01;   %low cut frequency in Hz 0.01 %play with this number matters a lot
Fc2 = 0.1;   %intend 0.1 Hz but consider transient band
N = order;
beta = 0.005;
win = kaiser(N+1, beta); %using kaiser window
flag = 'scale';  % Sampling Flag

% Calculate the coefficients using the FIR1 function.
b  = fir1(N, [Fc1, Fc2]/(Fs/2), 'bandpass', win, flag);

fvtool(b,1,'Fs',Fs)
Hd = dfilt.dffir(b);
cc_on = 1;% 0,w/o cc | 1,w/ cc
% apply filters
switch cc_on
    case 0
        cortical_depth_map = norm_cm;
    case 1
        cortical_depth_map = norm_cc;
end
for i = 1:size(fmri,3)
%     cortical_depth_map = norm_cm{i};
    filtered_cortical_depth_map = zeros(size(cortical_depth_map{i}));
    nVoxel = size(cortical_depth_map{i},1);
    nTime = size(filtered_cortical_depth_map,2);
    for voxel_num = 1 : nVoxel

        %prevent the signal fluctuation in front and behind part after filtering
        temp_filtered_cortical_depth_map1(voxel_num,:)=fliplr((cortical_depth_map{i}(voxel_num,2:N+1)));  % maintain continuity in level and slope
        temp_filtered_cortical_depth_map2(voxel_num,:)=fliplr((cortical_depth_map{i}(voxel_num,end-N:end-1)));

        filtered_cortical_depth_map_new(voxel_num,1:N)                 = temp_filtered_cortical_depth_map1(voxel_num,1:N);
        filtered_cortical_depth_map_new(voxel_num,N+1:nTime+N)         = (cortical_depth_map{i}(voxel_num,1:nTime));
        filtered_cortical_depth_map_new(voxel_num,nTime+N+1:nTime+2*N) = temp_filtered_cortical_depth_map2(voxel_num,1:N);

        %                         filtered_cortical_depth_map(voxel_num,:) =  filter(Hd,squeeze(cortical_depth_map(voxel_num,:)));
        filtered_cortical_depth_map_new(voxel_num,:) =  filter(Hd,squeeze(filtered_cortical_depth_map_new(voxel_num,:)));
        delay = mean(grpdelay(Hd));
        filtered_cortical_depth_map_new(voxel_num,:) = circshift(filtered_cortical_depth_map_new(voxel_num,:),(-1)*delay,2);
    end
    filtered_cortical_depth_map(:,1:nTime) =  filtered_cortical_depth_map_new(:,N+1:N+nTime);
    
    cortical_depth_map_all(:,:,i) = abs(filtered_cortical_depth_map);
end
% filtered fmri time course plot
for i = 1:size(fmri,3)
    slice_filter{i} = filter(Hd,fmri_tc{i});
    slice_mean_filter{i} = mean(abs(slice_filter{i}),1);
    slice_zscore_filter{i} = zscore(slice_mean_filter{i});
end
for i = 1:size(fmri,3)
    h = figure; 
%     set(h, 'WindowStyle', 'docked', 'NumberTitle', 'off', 'Name', ['filtered fMRI time course slice#', num2str(i)])
%     subplot(2,1,i);
    plot(1:6400,slice_zscore_filter{i},'b-');
    xlim([0,6400]);
    xlabel('time(sec)');
    ylabel('fMRI a.u.');
    title('averaged fmri filtered in 0.01-0.1Hz');
    set(gcf,'Position',[500 500 980 300 ]);
end
%% epoch-wise normalized colormap(cortical temporal BOLD map)
% cortical_depth_map = temp_cm_cut(:,:,3);


clear cortical_depth_map;
cortical_depth_map = cortical_depth_map_all;
cortical_depth_map
nY =32;
pre_stim = 10;
for i = 1:size(fmri,3)
    for nRepeat = 1 : nY %32

        SI_base_mean = mean(abs(cortical_depth_map(:,1+nSlice*(nRepeat-1):pre_stim+nSlice*(nRepeat-1),i)),2); % 1s off
        SI_i = abs(cortical_depth_map(:,1+nSlice*(nRepeat-1):nSlice*(nRepeat),i));% single epoch(1/32) response colormap
        for t = 1 : size(SI_i,2)
            SI_percentage(:,t) = (SI_i(:,t) - SI_base_mean(:,1))./(SI_base_mean(:,1));% percentage change within single epoch 
        end
        Threshold_SI_MAP(:,1+nSlice*(nRepeat-1):nSlice*(nRepeat)) = SI_percentage;
    end

    norm_cortical_map(:,:,i) = abs(Threshold_SI_MAP(:,:))./max(max(abs(Threshold_SI_MAP(:,:))));

    norm_cortical_map_temp(:,:,i) = abs(cortical_depth_map(:,:,i))./max(max(abs(cortical_depth_map(:,:,i))));
    Threshold_SI_MAP = abs(cortical_depth_map(:,:,i))./mean(norm_cortical_map_temp(:,:,i),2);
    norm_cortical_map(:,:,i) = abs(Threshold_SI_MAP(:,:))./max(max(abs(Threshold_SI_MAP(:,:))));
end
% plot colormap of BOLD 
TA = 640;
for i = 1:size(fmri,3)
    h = figure; set(h, 'WindowStyle', 'docked', 'NumberTitle', 'off', 'Name', ['fMRI time course slice', num2str(i)])
    imagesc(norm_cortical_map(:,:,i));
    cmax = max(max(norm_cortical_map(:,:,i)));
    cmin = min(min(norm_cortical_map(:,:,i)));
    colormap(jet);
    colorbar;
    xlim([0 TA*(1/TR)])
    xticks(0:TA/4*(1/TR):TA*(1/TR))
    xticklabels({'0s','160s','320s','480s','640s'});
    caxis([0,cmax])% i've no idea, but 0.3 makes colormap prettier :)
    set(gcf,'Position',[500 500 980 300 ]);
    saveas(h,'filtered timecourse vs depth.jpg');
end

%% Time course of ave fmri
for i = 1:size(fmri,3)
h = figure; 
% set(h, 'WindowStyle', 'docked', 'NumberTitle', 'off', 'Name', ['fMRI time course slice', num2str(i)])
mean_slice{i} = squeeze(mean(abs(fmri_tc{i}),1));
plot(1:6400,zscore(mean_slice{i}));
xlim([0,6400]);
xticks(0:6400/4*(1/TR):6400*(1/TR))
xlabel('time(s)');
xticklabels({'0 sec','160 sec','320 sec','480 sec','640 sec'})
ylabel('fMRI a.u.');
set(gcf,'Position',[500 500 980 300 ]);
end



% bpfit = designfilt('bandpassfir', 'FilterOrder', 20, 'CutoffFrequency1', 0.01, 'CutoffFrequency2', 0.1, 'SampleRate', 10);
for i = 1:2
figure;
[b,a] = butter(1,[0.01 5]/Fs/2,'bandpass');
temp_filt1 = filtfilt(b,a,mean(fmri_tc{i},1));
plot(1:6400,zscore(temp_filt1));
xlim([0,6400]);
ylim([-4,4]);
xticks(0:6400/4*(1/TR):6400*(1/TR))
xlabel('time(s)');
xticklabels({'0 sec','160 sec','320 sec','480 sec','640 sec'})
ylabel('fMRI a.u.');
set(gcf,'Position',[500 500 980 300 ]);
end
