%% Initializing parameters
addpath('C:\Users\wangqi\Documents\Lab\Demo\Fid2plot') % the path where Brucker functions were stored
path = 'D:\16072020_yi\'; % targeted path of data 
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
cc_on = 1;% 0,w/o cc | 1,w/ cc

% sampling rates;
TR = 0.1;
Fs = 1/TR;

% scans = [44,46,50];
scans = [91];

state = 'evoked'; % 'evoked' | 'rs'

% I/O path:

mkdir([path,'\figs']);
save_path=[path,'\figs'];

%% Read raw fid
for ir = 1:length(scans)
        path_fmri = [path,num2str(scans(ir)),'\fid'];
        [res, p, pathstr] = BrReadImage_AutoRun(path_fmri, params);
        res = squeeze(res);
        fmri(ir,:,:,:,:) = squeeze(res);
end
%% Average trails(evoked only)
% if state =='evoked'
% fmri = squeeze(mean(fmri,1));
% end
%% Formating fid
% fmri = squeeze(res);

nslice = size(fmri,4);
for ir = 1:length(scans)
    for is = 1:nslice
        fmri_slice{ir}{is} = squeeze(fmri(ir,:,:,is,:));
        slice_fftc{is} = fft2c(fmri_slice{ir}{is});
        slice_ifft{is} = ifft2c(slice_fftc{is});
        fmri_slice{ir}{is} = fftc(slice_ifft{is},2);
    %     [why not?] how does the `reshape` function organize matrix?
    %     fmri_slice{i} = reshape(slice_fftc{i},size(slice_fftc{i},1),(size(slice_fftc{i},2)*size(slice_fftc{i},3))); % [why?]by doing fft, the profile is more salient
        for k = 1:32
            for j = 1:200
                fmri_cm(:,k+32*(j-1)) = abs(fmri_slice{ir}{is}(:,k,j));
            end
        end
    fmri_tc{ir}{is} = fmri_cm; 
    
    end
end
[nX, nY, ndur] = size(squeeze(fmri(ir,:,:,is,:))); % nX = cortical depth | nY = epochs | nSlice = duration of epoch
t_fmri = 0:TR:(nY*ndur-1)*TR; % t fmri initialization
% tSNR plot
for ir = 1:length(scans)
    for is = 1:nslice
        h = figure; set(h, 'WindowStyle', 'docked', 'NumberTitle', 'off', 'Name', ['tSNR '])
        plot(linspace(1,nX,nX),abs(mean(fmri_tc{ir}{is},2)));
        xlabel('cortical depth(mm)');
        xlim([1,128]);
        ylabel('intensity');
        title(['Cortical surface trail#',num2str(scans(ir))]);
    end
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
cd(save_path)
% CC border still untargeted
for ir = 1:length(scans)
    for is = 1:nslice
        ave_tsnr(:,ir,is) = abs(mean(fmri_tc{ir}{is},2));
        ave_tsnr(:,ir,is) = ave_tsnr(:,ir,is) - min(ave_tsnr(:,ir,is));
        half_idx(:,ir,is) = find(ave_tsnr(:,ir,is)>(max(ave_tsnr(:,ir,is))/2),1);
%         half_idx(:,ir,is) = 44;
        half_idx_bot(:,ir,is) = half_idx(:,ir,is)+len_cortex; % bottom of cortex
        cc_idx_bot(:,ir,is) = half_idx_bot(:,ir,is)+depth_cc; % bottom of corpus callosum
    end
end

for ir = 1:length(scans)
    %
    for is = 1:nslice
        h = figure();
        subplot(1,2,1)
        imagesc(t_fmri,[],fmri_tc{ir}{is});
    %     imagesc(abs(fmri_slice{i}));
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
for ir = 1:length(scans)
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

for ir = 1:length(scans)
    for is = 1:nslice
        figure(),
        plot(1:size(tsnr_temp,1),squeeze(tsnr_temp(:,ir,is)));
        xlim([1,size(tsnr_temp,1)]);
        title('tSNR of cortex');
    end
end
%% cut off 
for ir = 1:length(scans)
    for is = 1:nslice
        norm_cm{ir}{is} = fmri_cm_cut{ir}{is} ./ max(max(fmri_cm_cut{ir}{is}));
        norm_cc{ir}{is} = fmri_cc_cut{ir}{is} ./ max(max(fmri_cc_cut{ir}{is})); % corpus callosum included
        figure(),
        cmin = min(min(norm_cm{ir}{is}));
        cmax = max(max(norm_cm{ir}{is}));
        imagesc(t_fmri,[],squeeze(norm_cm{ir}{is}));
        xlabel('time(s)');
        ylabel('cortical voxels');
        colormap(jet);
        caxis([cmin,cmax]);
        set(gcf,'Position',[500,500,600,300]);
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
        cc_on = 1;
        disp('filtering without Corpus callosum...')
    case 1
        cortical_depth_map = norm_cc;
        save_name_2 = 'temperal-spatial fmri map(with CC.)';
        cc_on = 0;
        disp('filtering with Corpus callosum...')
end 
% dimension of vairables, dependent of signal depth
filtered_cortical_depth_map = zeros(size(cortical_depth_map{ir}{is}));
nVoxel = size(cortical_depth_map{ir}{is},1); % num of voxels
nTime = size(filtered_cortical_depth_map,2); % num of time points
cortical_depth_map_all = zeros(nVoxel,nTime,size(fmri,3));
for ir = 1:length(scans)
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
for ir = 1:length(scans)
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

for ir = 1:length(scans)
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

        caxis([0,cmax])% i've no idea, but 0.3 makes colormap prettier :)
        set(gcf,'Position',[500 500 980 300 ]);
        saveas(h,[save_name_2,'trail# ',num2str(scans(ir)),' slice# ',num2str(is),'.jpg']);
    end
end
% clear SI_percentage Threshold_SI_MAP norm_cortical_map norm_cortical_map_temp
%% Time course of ave fmri
switch cc_on
    case 1 
        title_name = 'fmri time course including cc';
        title_name_filted = 'fmri time course including CC(0.01-1Hz)';
        fmri_timecourse = fmri_cc_cut;
        cc_on = 0;
    case 0
        title_name = 'fmri time course';
        title_name_filted = 'fmri time course(0.01-0.1Hz)';
        fmri_timecourse = fmri_cm_cut;
        cc_on = 1;
end
for ir = 1:length(scans)
    for is = 1:nslice
    h = figure; 
    % set(h, 'WindowStyle', 'docked', 'NumberTitle', 'off', 'Name', ['fMRI time course slice', num2str(i)])
    mean_slice_1(ir,is,:) = squeeze(mean(abs(fmri_timecourse{ir}{is}),1));
    plot(t_fmri,zscore(squeeze(mean_slice_1(ir,is,:))));
    xlim([0,t_fmri(end)]);
    % xticks(0:6400/4*(1/TR):6400*(1/TR));
    xlabel('time(s)');
    % xticklabels({'0 sec','160 sec','320 sec','480 sec','640 sec'})
    ylabel('fMRI a.u.');
    title([title_name,'trail#',num2str(scans(ir)),' slice#',num2str(is)]);
    set(gcf,'Position',[500 500 980 300 ]);
    box off;
    saveas(h,[title_name,'trail#',num2str(scans(ir)),' slice#',num2str(is),'.jpg']);
    end
end


% filtering fmri
bpf_low = 0.01;
bpf_high = 0.1;
[b,a] = butter(1,[bpf_low bpf_high]/(Fs/2),'bandpass');
len = size(fmri_timecourse{ir}{is},2);
for ir = 1:length(scans)
    for is = 1:nslice
        fmri_ave{ir}{is} = mean(fmri_timecourse{ir}{is},1);
        temp = [fliplr(fmri_ave{ir}{is}(:,1:len/4)),fmri_ave{ir}{is},fliplr(fmri_ave{ir}{is}(:,end-len/4+1:end))];
        temp_filt1 = filtfilt(b,a,temp);
        fmri_ave_f{ir}{is} = temp_filt1(:,len/4+1:end-len/4);
    end
end
for ir = 1:length(scans)
    for is = 1:nslice
        h = figure;
        plot(t_fmri,zscore(fmri_ave_f{ir}{is}));
        xlim([0,t_fmri(end)]);
        ylim([-4,4]);
        xlabel('time(s)');
        ylabel('fMRI a.u.');
        grid on;
        box off;
        title([title_name_filted,'trail#',num2str(scans(ir)),' slice#',num2str(is)]);
        set(gcf,'Position',[500 500 980 300 ]);
        saveas(h,[title_name_filted,'trail#',num2str(scans(ir)),' slice#',num2str(is),'.jpg']);
    end
end
% spectral analysis
for is = 1:2
    test(is,:) = squeeze(mean_slice_1(ir,is,:));
    test(is,:) = test(is,:) - mean(test(is,:));
    h = figure;
    [pxx,f] = pwelch(test(is,:),size(test,2),0,[],1/TR);
    plot(f(f>0&f<0.12),pxx(f>0&f<0.12,:));
    xlim([f(2),0.12]);
    ylim([0,50]);
    xlabel('frequency(Hz)');
    ylabel('intensity dB');
    title(['spectrum of slice#',num2str(is)]);
    saveas(h,['spectrum of slice#',num2str(is)]);
end
%% Power spectrogram[NOT YET COMPLETED]
% NFFT = 2^12;
% window = 1*Fs; % 1sec of windowing
% overlap = 0;
% for i = 1:size(fmri,3)
%  h = figure;
% [S,F,T] = spectrogram(mean_slice_1(i,:),window,overlap,NFFT,Fs,'yaxis');
% f_sample = [0:(Fs/2)/NFFT:(length(F)-1)*(Fs/2)/NFFT];
% imagesc(T,f_sample(f_sample<1),abs(S));
% set(gca,'YDir','normal');
% xlabel('time(s)');c
% ylabel('Frequency(Hz)');
% title('fMRI spectrogram(win = 1 sec, no overlap)');
% box off;
% end