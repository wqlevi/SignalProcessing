% IF USING THIS CODE ON WINDOWS CHANGE ALL '/' to '\' in paths

clear

path_data = 'C:\Users\wangqi\Documents\Lab\Data\05302019_Hang\05302019_analyzed';

mode_paradigm = 'task'; % 'rs' | 'task'
mode_cortex = 'voxels'; % 'voxels' | 'layers'

do_convolve = 0; % do you want to convolve all ca feature signals with an HRF kernel?

%%
TR50 = 0; % the TR is either 0.05 or 0.1 s; set TR50 to 1 if its 0.05

chan_ca = 7; % channel with calcium data

ws = 2^14; % spectrogram window length. Needs to be specified beforehand because of specific matching procedure.

n_plots = 21;

axis_width = 2;
font_size=10;

%% scan params
if TR50 == 1
    TR = 0.05; %s
    b4trig = 20; %prestimbaseline (was 20 before SChoi preprocessing)
    total_scan_time = 640; %second
else
    TR = 0.1; %s
    b4trig = 10; %prestimbaseline (was 10 before SChoi preprocessing)
    total_scan_time = 640; %second
end
n_slices = 3;
epoch_len = 20; % seconds
n_epochs = total_scan_time / epoch_len;
fmri_dummy = ones(total_scan_time*(1/TR), 1);

%% runs list
switch mode_paradigm
    case 'task'              
%           runs = [73 75 77 79 81]; % task-fMRI
        runs = [16 17 20 23 25 27 29 31 33 35 39 41 43 45 47 49 51 53 55 59 61 63 65 67 69 71 73 75 77 79 81];
%         runs = [15];
        path_mode = 'C:\Users\wangqi\Documents\Lab\Data\05302019_Hang\05302019_ca';
%         path_suffix = '_TK.mat';
        path_suffix = '.mat';
    case 'rs'
%         runs = [17];
%           runs = [17 19  21  22 24 26 27 32 33 36 38 40 42 44 48 50 51 53 55 57]; %10062018, rs-fMRI
%         runs = [51]; % rs-fMRI
        runs = [15 19 22 24 26 28 30 32 34 36 37 38 40 42 44 46 48 50 52 54 56 57 58 60 62 64 66 68 70 72 74 76 78 80];
        path_mode = 'C:\Users\wangqi\Documents\Lab\Data\05302019_Hang\05302019_ca';
%         path_suffix = '_RS.mat';
        path_suffix = '.mat';
end

if length(runs)<n_plots
    n_plots = length(runs);
end

%% load signals
disp('load signals')
for ir = 1 : length(runs)
    %% load 3 slice data for all runs
    for is = 1:n_slices
%         tmp = load([path_data, '1.PI1/', num2str(runs(ir)), '/Results_Slice', num2str(is), ...
%             '/line_scanning_data_s', num2str(is), '.mat']);
   tmp = load([path_data, '\', num2str(runs(ir)), '\Results_Slice', num2str(is), ...
            '\line_scanning_data_s', num2str(is), '.mat']);
        fmri{is}{ir} = abs(tmp.total_cortical_depth_map)';
    end
    %% get calcium signal exactly matching the fmri signal
%     data_orig = load([path_data, path_mode, 'scan_',num2str(runs(ir)), path_suffix]);
    path_data_ca='C:\Users\wangqi\Documents\Lab\Data\05302019_Hang\05302019_ca';
    data_orig = load([path_data_ca, '\scan_',num2str(runs(ir)), path_suffix]);
    [data_match, fmri_dummy, beg, fin] = match_acq_fmri(data_orig, fmri_dummy, TR, b4trig);
    ca_match{ir} = -(data_match.channels{chan_ca}.data)'; % Reversed raw calcium data after match
    fs = (data_match.channels{chan_ca}.samples_per_second); % Freq. of calcium sampling
    
    %% get calcium slightly longer than BOLD so that spectrogram times are aligned to fMRI times
    if beg(chan_ca) > ws/2 - TR*fs/2
        beg_long = beg(chan_ca) - ws/2 + TR*fs/2;
    else
        beg_long = 1;
        disp(['Run ',num2str(run),' - signal too short at the beginning. Alignment wont be perfect.'])
    end
    if fin(chan_ca) + ws/2 - TR*fs/2 < length(data_orig.channels{chan_ca}.data)
        fin_long = fin(chan_ca) + ws/2 - TR*fs/2;
    else
        fin_long = length(data_orig.channels{chan_ca}.data);
        disp(['Run ',num2str(run),' - signal too short at the end. Alignment wont be perfect.'])
    end
    ca_long{ir} = -data_orig.channels{chan_ca}.data(beg_long:fin_long)';
end

%% demean Ca
for ir = 1 : length(runs)
    ca_match{ir} = ca_match{ir} - mean(ca_match{ir});
    ca_long{ir} = ca_long{ir} - mean(ca_long{ir});
    % DON'T DEMEAN FMRI HERE (NOT BEFORE FINDING CORTEX BORDERS)
    % %     fmri_s1{ir} = bsxfun(@minus, fmri_s1{ir}, mean(fmri_s1{ir}, 1));
    % %     fmri_s2{ir} = bsxfun(@minus, fmri_s2{ir}, mean(fmri_s2{ir}, 1));
    % %     fmri_s3{ir} = bsxfun(@minus, fmri_s3{ir}, mean(fmri_s3{ir}, 1));
end

%% generate time
t_fmri = [0 : TR : (length(fmri_dummy)-1)*TR]';
t_ca = [0 : 1/fs : (length(ca_match{1})-1)/fs]';

%% find the corex boundary
disp('find cortex boundary')

for is = 1:n_slices
    fmri_ccat{is} = cell2mat(fmri{is}');
    mean_slice = mean(fmri_ccat{is}, 1);
    mean_slice = mean_slice - min(mean_slice);
    half_intens = max(mean_slice)/2;
    i_half_ccat(is) = find(mean_slice>half_intens, 1);
end

for ir = 1:length(runs)
    for is = 1:n_slices
        mean_slice = mean(fmri{is}{ir}, 1);
        mean_slice = mean_slice - min(mean_slice);
        
        half_intens = max(mean_slice)/2;
        i_half(is,ir) = find(mean_slice>half_intens, 1);
        
        % shift the cortex lower by one voxel
        i_half(is,ir) = i_half(is,ir) + 1;
    end
end

%% plot cortex borders
for is =1:n_slices
    h = figure; set(h, 'WindowStyle', 'docked', 'NumberTitle', 'off', 'Name', ['ctx border ', num2str(is)])
    imagesc(1:size(fmri_ccat{is}, 1), 1:size(fmri_ccat{is}, 2), fmri_ccat{is}'), hold on
    plot(1:size(fmri_ccat{is}, 1), ones(1, size(fmri_ccat{is}, 1))*i_half_ccat(is), 'w', 'LineWidth', 2)
    plot(1:size(fmri_ccat{is}, 1), ones(1, size(fmri_ccat{is}, 1))*(i_half_ccat(is)+41), 'w', 'LineWidth', 2)
    for ir = 1:length(runs)
        len = size(fmri{is}{ir}, 1);
        plot(1+(ir-1)*len:ir*len, ones(1, len)*i_half(is,ir), 'r', 'LineWidth', 2)
        plot(1+(ir-1)*len:ir*len, ones(1, len)*(i_half(is,ir)+41), 'r', 'LineWidth', 2)
    end
    box off
    title(['cortex borders, all trials, slice ', num2str(is)])
end

%% extract cortical layer/voxel signals
% Spatial resolution, 0.05 mm
for is = 1:n_slices
    for ir = 1:length(runs)
        ctx_sig_tmp = fmri{is}{ir}(:, i_half(is,ir)+1 : i_half(is,ir)+40);
        switch mode_cortex
            case 'layers'
                ctx{is}{ir}(:, 1) = mean(ctx_sig_tmp(:, 1:3), 2);    % L1, 0~0.15 mm, averaged voxel
                ctx{is}{ir}(:, 2) = mean(ctx_sig_tmp(:, 4:11), 2);   % L2/3, 0.15~0.55 mm, averaged voxel
                ctx{is}{ir}(:, 3) = mean(ctx_sig_tmp(:, 12:16), 2);  % L4, 0.55~0.8 mm, averaged voxel
                ctx{is}{ir}(:, 4) = mean(ctx_sig_tmp(:, 17:27), 2);  % L5, 0.8~1.4 mm, averaged voxel
                ctx{is}{ir}(:, 5) = mean(ctx_sig_tmp(:, 28:40), 2);  % L6, 1.4~2.0 mm, averaged voxel
            case 'voxels'
                ctx{is}{ir} = ctx_sig_tmp;
        end
    end
end

%% deman fMRI
% % for is=1:n_slices
% %     for ir=1:length(runs)
% %        ctx{is}{ir} = ctx{is}{ir} - mean(ctx{is}{ir}, 2);
% %     end
% % end

%% alternatively variance normalize and demean
% % for is=1:n_slices
% %     for ir=1:length(runs)
% %        ctx{is}{ir} = zscore(ctx{is}{ir}, [], 1);
% %     end
% % end

%% regress polynomial of degree n from fMRI data
for is=1:n_slices
    for ir=1:length(runs)
        for iv = 1:size(ctx{is}{ir}, 2)
            p = polyfit(t_fmri, ctx{is}{ir}(:, iv), 3);
            ctx{is}{ir}(:, iv) = ctx{is}{ir}(:, iv) - polyval(p, t_fmri);
        end
    end
end

%% REMOVE POLY CA

%% filter fMRI
if strcmp(mode_paradigm, 'rs')
    disp('filter fMRI')
    [b,a] = butter(1, [0.01, 0.1]/((1/TR)/2), 'bandpass');
    for is = 1:n_slices
        for ir = 1:length(runs)
            ctx{is}{ir} = filtfilt(b,a,ctx{is}{ir});
        end
    end
end

%% fMRI - between slice lags
for ir = 1:length(runs)
    [xxx12(ir, :), lag] = xcorr(zscore(mean(ctx{1}{ir}, 2)), zscore(mean(ctx{2}{ir}, 2)), 100, 'coeff');
    [xxx13(ir, :), lag] = xcorr(zscore(mean(ctx{1}{ir}, 2)), zscore(mean(ctx{3}{ir}, 2)), 100, 'coeff');
    [xxx23(ir, :), lag] = xcorr(zscore(mean(ctx{2}{ir}, 2)), zscore(mean(ctx{3}{ir}, 2)), 100, 'coeff');
end
h = figure; set(h, 'WindowStyle', 'docked', 'NumberTitle', 'off', 'Name', ['fmri xcorr s12'])
plot(lag*TR, mean(xxx12, 1))
h = figure; set(h, 'WindowStyle', 'docked', 'NumberTitle', 'off', 'Name', ['fmri xcorr s13'])
plot(lag*TR, mean(xxx13, 1))
h = figure; set(h, 'WindowStyle', 'docked', 'NumberTitle', 'off', 'Name', ['fmri xcorr s23'])
plot(lag*TR, mean(xxx23, 1))

%% compute spectrograms & extract freq. band signals
disp('compute spectrograms')
for ir = 1:length(runs)
    win = hanning(ws);
    ovrl = ws - (length(ca_long{ir})-ws)/(length(fmri_dummy)-1);
    nfft = 2^15;
    [S,F,T] = spectrogram(ca_long{ir}, win, ovrl, nfft, fs);
    
    if ir <= n_plots
        f_int = F < 15 & F>0.5;
        h = figure; set(h, 'WindowStyle', 'docked', 'NumberTitle', 'off', 'Name', ['spectrogram r', num2str(runs(ir))])
        %         imagesc(T, F(f_int), log10(abs(S(f_int, :)))), box off
        imagesc(T, F(f_int), abs(S(f_int, :))), box off
        ax = gca;
        ax.YDir = 'normal';
        title(['spectrogram, trial ', num2str(runs(ir))]);
        xlabel('Time(s)');ylabel('Frequence(Hz)');colorbar;
        box off;
    end
    
    f_low = 1; f_high = 5;
    ca_freq_15{ir} = mean(abs(S(F<=f_high&F>=f_low, :)))';
    %     f_high = 1;
    %     ca_freq_01{ir} = mean(abs(S(F<=f_high, :)))';
    %     f_high = 5;
    %     ca_freq_05{ir} = mean(abs(S(F<=f_high, :)))';
end

%% filter freq. bands
disp('filter freq')
f_low = 0.01;
f_high_1 = 0.1;
[b1, a1] = butter(1, [f_low f_high_1]/((1/TR)/2), 'bandpass');
for ir = 1:length(runs)
    len = length(ca_freq_15{ir});
    freq_band_tmp = [fliplr(ca_freq_15{ir}(1:round(len/4))); ca_freq_15{ir}; ...
        fliplr(ca_freq_15{ir}(end-round(len/4)+1:end))];
    ca_freq_15{ir} = filtfilt(b1, a1, freq_band_tmp);
    
    %     freq_band_tmp = [fliplr(ca_freq_01{ir}(1:round(len/4))); ca_freq_01{ir}; ...
    %         fliplr(ca_freq_01{ir}(end-round(len/4)+1:end))];
    %     ca_freq_01{ir} = filtfilt(b1, a1, freq_band_tmp);
    %
    %     freq_band_tmp = [fliplr(ca_freq_05{ir}(1:round(len/4))); ca_freq_05{ir}; ...
    %         fliplr(ca_freq_05{ir}(end-round(len/4)+1:end))];
    %     ca_freq_05{ir} = filtfilt(b1, a1, freq_band_tmp);
    
    ca_freq_15{ir} = ca_freq_15{ir}(round(len/4)+1 : end-round(len/4));
    %     ca_freq_01{ir} = ca_freq_01{ir}(round(len/4)+1 : end-round(len/4));
    %     ca_freq_05{ir} = ca_freq_05{ir}(round(len/4)+1 : end-round(len/4));
end

%% filter ca
disp('filter Ca')
f_low = 0.01;
f_high_1 = 0.1;
f_high_2 = 20;
% f_high_3 = 0.03;

[b1, a1] = butter(1, [f_low f_high_1]/(fs/2), 'bandpass');
[b2, a2] = butter(1, [f_low f_high_2]/(fs/2), 'bandpass');
% [b3, a3] = butter(1, f_high_3/(fs/2), 'low');

for ir = 1:length(runs)
    len = length(ca_match{ir});
    ca_tmp = [fliplr(ca_match{ir}(1:round(len/4))); ca_match{ir}; ...
        fliplr(ca_match{ir}(end-round(len/4)+1:end))];
    ca_f1{ir} = filtfilt(b1, a1, ca_tmp);
    ca_f2{ir} = filtfilt(b2, a2, ca_tmp);
    %     ca_f3{ir} = filtfilt(b3, a3, ca_tmp);
    %
    ca_f1{ir} = ca_f1{ir}(round(len/4)+1 : end-round(len/4));
    ca_f2{ir} = ca_f2{ir}(round(len/4)+1 : end-round(len/4));
    %     ca_f3{ir} = ca_f3{ir}(round(len/4)+1 : end-round(len/4));
    %
    ca_f1_ds{ir} = downsample(ca_f1{ir}, fs*TR);
    ca_f2_ds{ir} = downsample(ca_f2{ir}, fs*TR);
    %     ca_f3_ds{ir} = downsample(ca_f3{ir}, fs*TR);
end

%% convolve all ca sigs
% based on Patricia's script
if do_convolve
    disp('convolve Ca sigs')
    
    % build HRF
    time=1:1:30; % HRF duration
%     T0=0; n=3; lamda=2;
    T0=0; n=5; lamda=1;
    hrf=((time-T0).^(n-1)).*exp(-(time-T0)/lamda)/((lamda^n)*factorial(n-1));
    h=hrf/max(hrf); %figure,plot(h)
    
    % resample gamma/HRF to TR
    oldsr=1; % from seconds..
    newsr=1/TR; % ..to 1/TR
    [P,Q] = rat(newsr/oldsr);
    rs_h=resample(h,P,Q); % figure,plot(rs_h);
    
    % convolve
    for ir =1:length(runs)
        ca_f1_ds{ir} = conv(rs_h, ca_f1_ds{ir});
        ca_f1_ds{ir} = ca_f1_ds{ir}(1:length(t_fmri));
        
        ca_f2_ds{ir} = conv(rs_h, ca_f2_ds{ir});
        ca_f2_ds{ir} = ca_f2_ds{ir}(1:length(t_fmri));
        
        ca_freq_15{ir} = conv(rs_h, ca_freq_15{ir});
        ca_freq_15{ir} = ca_freq_15{ir}(1:length(t_fmri));
    end
end

%% plot HRF convolution kernel
% h = figure; set(h, 'WindowStyle', 'docked', 'NumberTitle', 'off', 'Name', ['HRF kernel'])
% plot(rs_h)

%% x-corr - filt
clear xc_f1 xc_f2 xc_f3 xc_15 xc_05 xc_f01
disp('x-corr ca vs. fmri')

switch mode_paradigm
    case 'task'
        xc_lag = 100;
    case 'rs'
        xc_lag = 300;
end

for is = 1:n_slices
    for ir = 1:length(runs)
        for iL = 1: size(ctx{is}{ir}, 2)
            [xc_f1{is}(ir, :, iL), lag] = xcorr(zscore(ca_f1_ds{ir}), zscore(ctx{is}{ir}(:, iL)), xc_lag, 'coeff');
            [xc_f2{is}(ir, :, iL), lag] = xcorr(zscore(ca_f2_ds{ir}), zscore(ctx{is}{ir}(:, iL)), xc_lag, 'coeff');
            %             [xc_f3{is}(ir, :, iL), lag] = xcorr(zscore(ca_f3_ds{ir}), zscore(ctx{is}{ir}(:, iL)), xc_lag, 'coeff');
            
            [xc_15{is}(ir, :, iL), lag] = xcorr(zscore(ca_freq_15{ir}), zscore(ctx{is}{ir}(:, iL)), xc_lag, 'coeff');
            %             [xc_05{is}(ir, :, iL), lag] = xcorr(zscore(ca_freq_05{ir}), zscore(ctx{is}{ir}(:, iL)), xc_lag, 'coeff');
            %             [xc_01{is}(ir, :, iL), lag] = xcorr(zscore(ca_freq_01{ir}), zscore(ctx{is}{ir}(:, iL)), xc_lag, 'coeff');
        end
    end
    
    % cc at 0 lag
    cc_f1{is} = squeeze(xc_f1{is}(:, lag==0, :));
    cc_f2{is} = squeeze(xc_f2{is}(:, lag==0, :));
    %     cc_f3{is} = squeeze(xc_f3{is}(:, lag==0, :));
    
    cc_15{is} = squeeze(xc_15{is}(:, lag==0, :));
    %     cc_05{is} = squeeze(xc_05{is}(:, lag==0, :));
    %     cc_01{is} = squeeze(xc_01{is}(:, lag==0, :));
    
    %
    [~, im_f1{is}] = max(xc_f1{is}, [], 2);
    [~, im_f2{is}] = max(xc_f2{is}, [], 2);
    %     [~, im_f3{is}] = max(xc_f3{is}, [], 2);
    
    [~, im_15{is}] = max(xc_15{is}, [], 2);
    %     [~, im_05{is}] = max(xc_05{is}, [], 2);
    %     [~, im_01{is}] = max(xc_01{is}, [], 2);
end

%% plot ca filtered 0.01-20 Hz
for ir = 1:length(runs)
    h = figure; set(h, 'WindowStyle', 'docked', 'NumberTitle', 'off', 'Name', ['ca 0.01-20 r', num2str(runs(ir))])
    plot(t_ca, ca_f2{ir})
    box off, axis tight
end

%% lag times
for ir = 1:length(runs)
    h = figure; set(h, 'WindowStyle', 'docked', 'NumberTitle', 'off', 'Name', ['lag times r', num2str(runs(ir))])
    imagesc(1:n_slices, 1:length(im_f1{1}), squeeze([im_f1{1}(ir,:,:); im_f1{2}(ir,:,:); im_f1{3}(ir,:,:)])'), colorbar, colormap(jet)
    %     imagesc(1, 1:length(im), squeeze(im)'), colorbar, colormap(jet)
end

%% corr across voxels and slices
for ir = 1:length(runs)
    h = figure; set(h, 'WindowStyle', 'docked', 'NumberTitle', 'off', 'Name', ['cc f1 r', num2str(runs(ir))])
    imagesc(1:n_slices, 1:length(cc_f1{1}), [cc_f1{1}(ir,:,:); cc_f1{2}(ir,:,:); cc_f1{3}(ir,:,:)]'), colorbar, colormap(jet)
    %     imagesc(1, 1:length(im), squeeze(im)'), colorbar, colormap(jet)
end
cx = [mean(cc_f1{1}, 1); mean(cc_f1{2}, 1); mean(cc_f1{3}, 1)]';
h = figure; set(h, 'WindowStyle', 'docked', 'NumberTitle', 'off', 'Name', ['cc f1 mean'])
imagesc(1:n_slices, 1:length(cc_f1{1}), cx), colorbar, colormap(jet)
% % for ir = 1:length(runs)
% %     h = figure; set(h, 'WindowStyle', 'docked', 'NumberTitle', 'off', 'Name', ['cc r', num2str(runs(ir))])
% %     imagesc(1:n_slices, 1:length(cc_f2{1}), [cc_f2{1}(ir,:,:); cc_f2{2}(ir,:,:); cc_f2{3}(ir,:,:)]'), colorbar, colormap(jet)
% % %     imagesc(1, 1:length(im), squeeze(im)'), colorbar, colormap(jet)
% % end
for ir = 1:length(runs)
    h = figure; set(h, 'WindowStyle', 'docked', 'NumberTitle', 'off', 'Name', ['cc 15 r', num2str(runs(ir))])
    imagesc(1:n_slices, 1:length(cc_15{1}), [cc_15{1}(ir,:,:); cc_15{2}(ir,:,:); cc_15{3}(ir,:,:)]'), colorbar, colormap(jet)
    %     imagesc(1, 1:length(im), squeeze(im)'), colorbar, colormap(jet)
end
cx = [mean(cc_15{1}, 1); mean(cc_15{2}, 1); mean(cc_15{3}, 1)]';
h = figure; set(h, 'WindowStyle', 'docked', 'NumberTitle', 'off', 'Name', ['cc 15 mean'])
imagesc(1:n_slices, 1:length(cc_15{1}), cx), colorbar, colormap(jet)
% % cx = [mean(cc_05{1}, 1); mean(cc_05{2}, 1); mean(cc_05{3}, 1)]';
% % h = figure; set(h, 'WindowStyle', 'docked', 'NumberTitle', 'off', 'Name', ['cc r', num2str(runs(ir))])
% %     imagesc(1:n_slices, 1:length(cc_05{1}), cx), colorbar, colormap(jet)

%% lag curves at each voxel 2D
% % for is=1:n_slices
% %     [~, im] = max(xc_f1{is}, [], 2);
% %     lag_f1{is} = mean(squeeze(lag(im)), 1)*TR;
% % end
for ir = 1:length(runs)
    h = figure; set(h, 'WindowStyle', 'docked', 'NumberTitle', 'off', 'Name', ['lag 2D f1 r', num2str(runs(ir))])
    for is=1:n_slices
        subplot(1,3,is)
        imagesc(lag*TR, 1:size(xc_f1{is}, 3), squeeze(xc_f1{is}(ir,:,:))')
        cmin = min(min(squeeze(xc_f1{is}(ir,:,:))'));
        cmax = max(max(squeeze(xc_f1{is}(ir,:,:))'));
        title (['slice',num2str(is)]);
        xlabel('Times(s)');
        ylabel ('Voxels');
    end
    for is=1:n_slices
        subplot(1,3,is)
        colorbar
        caxis([cmin cmax])
    end
    
    h = figure; set(h, 'WindowStyle', 'docked', 'NumberTitle', 'off', 'Name', ['lag 2D f2 r', num2str(runs(ir))])
    for is=1:n_slices
        subplot(1,3,is)
        imagesc(lag*TR, 1:size(xc_f2{is}, 3), squeeze(xc_f2{is}(ir,:,:))')
        cmin = min(min(squeeze(xc_f2{is}(ir,:,:))'));
        cmax = max(max(squeeze(xc_f2{is}(ir,:,:))'));
        title (['slice',num2str(is)]);
        xlabel('Times(s)');
        ylabel ('Voxels');
    end
    for is=1:n_slices
        subplot(1,3,is)
        colorbar
        caxis([cmin cmax])
    end
    
    h = figure; set(h, 'WindowStyle', 'docked', 'NumberTitle', 'off', 'Name', ['lag 2D 15 r', num2str(runs(ir))])
    for is=1:n_slices
        subplot(1,3,is)
        imagesc(lag*TR, 1:size(xc_15{is}, 3), squeeze(xc_15{is}(ir,:,:))')
        cmin = min(min(squeeze(xc_15{is}(ir,:,:))'));
        cmax = max(max(squeeze(xc_15{is}(ir,:,:))'));
        title (['slice',num2str(is)]);
        xlabel('Times(s)');
        ylabel ('Voxels');
    end
    for is=1:n_slices
        subplot(1,3,is)
        colorbar
        caxis([cmin cmax])
    end
end

%% lag curves at each voxel 1D
% % for is=1:n_slices
% %     [~, im] = max(xc_f1{is}, [], 2);
% %     lag_f1{is} = mean(squeeze(lag(im)), 1)*TR;
% % end
for ir = 1:length(runs)
    h = figure; set(h, 'WindowStyle', 'docked', 'NumberTitle', 'off', 'Name', ['lag curves f1 r', num2str(runs(ir))])
    for is=1:n_slices
        subplot(1,3,is)
        plot(lag*TR, squeeze(xc_f1{is}(ir,:,:))')
        title (['slice',num2str(is)]);
        xlabel('Times(s)');
        ylabel ('Cross correlation coefficient');
    end
    
    h = figure; set(h, 'WindowStyle', 'docked', 'NumberTitle', 'off', 'Name', ['lag curves f2 r', num2str(runs(ir))])
    for is=1:n_slices
        subplot(1,3,is)
        plot(lag*TR, squeeze(xc_f2{is}(ir,:,:))')
        title (['slice',num2str(is)]);
        xlabel('Times(s)');
        ylabel ('Cross correlation coefficient');

    end
   
    h = figure; set(h, 'WindowStyle', 'docked', 'NumberTitle', 'off', 'Name', ['lag curves 15 r', num2str(runs(ir))])
    for is=1:n_slices
        subplot(1,3,is)
        plot(lag*TR, squeeze(xc_15{is}(ir,:,:))')
        title (['slice',num2str(is)]);
        xlabel('Times(s)');
        ylabel ('Cross correlation coefficient');
        
% %       %Starting Qi's editing
        x_lag_time = lag*TR;
        y_coefficient = squeeze(xc_15{is}(ir,:,:))';
%         for index_voxels = 1:size(y_coefficient,1)
%             hold on;
%             grid on;
%             plot(x_lag_time,y_coefficient(index_voxels,:));
%             [pks(index_voxels,:),locs(index_voxels,:)]=findpeaks(y_coefficient(index_voxels,:),x_lag_time,'NPeaks',1);
%             plot(locs(index_voxels,:),pks(index_voxels,:),'x');
%             %             findpeaks(y_coefficient(index_voxels,:),x_lag_time);
%             sprintf('Voxels #%d',index_voxels)
%             sprintf('scans #%d',runs)
%             hold off;
%             name_xlxs = sprintf('Scan #%d.xlsx',runs);
%         end  
% %          xlswrite(name_xlxs,locs{is},is);% Excel:sheet n contains lag time of slice n
            
% doesn't wrok well
            x_lag_time = lag*TR;
            y_coefficient = squeeze(xc_15{is}(ir,:,:))';
            for index_voxels = 1:size(y_coefficient,1)
                hold on;
                grid on;
                [pks(index_voxels,ir,is),locs(index_voxels,ir,is)] = max(y_coefficient(index_voxels,:),[],2);
                locs(index_voxels,ir,is) = (locs(index_voxels,ir,is)-100)*TR;
                plot(x_lag_time,y_coefficient);
                plot(locs(index_voxels,ir,is),pks(index_voxels,ir,is),'x');
                hold off;               
            end
%             name_xlxs = sprintf('Scan #%d.xlsx',runs);
%             xlswrite(name_xlxs,locs(:,ir,is));
            

%%%%%%%%%fMRI ave%%%%%%%%%%%%%%%%%%%%
 avg = zeros(200,1);
 figure,
 for index_slice = 1:3

for i=1:32
avg = avg+mean(ctx{index_slice}{31}(1+(i-1)*200:i*200,:),2);
end
subplot(3,1,index_slice);
plot(avg);
xlabel('time(s)');
ylabel('amplitude');
grid on;
name_slice = sprintf('Averaged fmri signal in single scan(slice %d)',index_slice);
title(name_slice);
 end
 %%%%%%%%%%%%%ca_ave
 avg = zeros(200,1);
 f_ca_ave = figure(),
for i=1:32
% avg = avg+mean(ctx{index_slice}{31}(1+(i-1)*200:i*200,:),2);

avg = avg+mean(ca_freq_15{ir}(1+(i-1)*200:i*200,:),2);

end
plot(avg);
xlabel('time(s)');
ylabel('amplitude');
grid on;
% name_slice = sprintf('Averaged fmri signal in single scan(slice %d)',index_slice);
name_slice = sprintf('Averaged Calcium signal in single scan(%d)',ir);
title(name_slice);
saveas(f_ca_ave,name_slice,'png');
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       %Finished Qi's editing
    locs_cell = num2cell(locs,[1 2]);
    writematrix(locs_cell{is},'xcorr_Scan_TK_slice.xlsx','Sheet',is);%sheet for slice, x-axis for trails, y-axis for voxels
    end

end

%% mean voxel & ca response across all trials & epochs. For evoked only
if strcmp(mode_paradigm, 'task')
    epoch_len = length(t_fmri)/n_epochs;
    
    % fmri
    for is = 1:n_slices
        ctx_resp{is} = zeros(epoch_len, size(ctx{1}{1}, 2));
        for ir = 1:length(runs)
            for ie = 1:n_epochs
                ctx_resp{is} = ctx_resp{is} + ctx{is}{ir}(1+(ie-1)*epoch_len:ie*epoch_len, :);
            end
            ctx_resp{is} = ctx_resp{is} / (n_epochs*length(runs));
        end
    end
    
    % ca
    ca_f1_resp = zeros(epoch_len, 1);
    ca_f2_resp = zeros(epoch_len, 1);
    ca_15_resp = zeros(epoch_len, 1);
    for ir = 1:length(runs)
        for ie = 1:n_epochs
            ca_f1_resp = ca_f1_resp + ca_f1_ds{ir}(1+(ie-1)*epoch_len:ie*epoch_len, :);
            ca_f2_resp = ca_f2_resp + ca_f2_ds{ir}(1+(ie-1)*epoch_len:ie*epoch_len, :);
            ca_15_resp = ca_15_resp + ca_freq_15{ir}(1+(ie-1)*epoch_len:ie*epoch_len, :);
        end
    end
    ca_f1_resp = ca_f1_resp / (n_epochs*length(runs));
    ca_f2_resp = ca_f2_resp / (n_epochs*length(runs));
    ca_15_resp = ca_15_resp / (n_epochs*length(runs));
    
end

%% x-corr of mean responses
if strcmp(mode_paradigm, 'rs')
    for is = 1:n_slices
        for iL = 1: size(ctx{is}{ir}, 2)
            [xc_f1_mean{is}(:, iL), lag] = xcorr(zscore(ca_f1_resp), zscore(ctx_resp{is}(:, iL)), xc_lag, 'coeff');
            [xc_f2_mean{is}(:, iL), lag] = xcorr(zscore(ca_f2_resp), zscore(ctx_resp{is}(:, iL)), xc_lag, 'coeff');
            [xc_15_mean{is}(:, iL), lag] = xcorr(zscore(ca_15_resp), zscore(ctx_resp{is}(:, iL)), xc_lag, 'coeff');
        end
    end
end

%% plot mean voxel & fMRI curves
if strcmp(mode_paradigm, 'rs')
    h = figure; set(h, 'WindowStyle', 'docked', 'NumberTitle', 'off', 'Name', ['ctx resp meanEpoch'])
    for is = 1:3
    subplot(3,1,is)
    plot(ctx_resp{is})
    end

    h = figure; set(h, 'WindowStyle', 'docked', 'NumberTitle', 'off', 'Name', ['ca resp meanEpoch'])
    subplot(3,1,1)
    plot(ca_f1_resp)
    subplot(3,1,2)
    plot(ca_f2_resp)
    subplot(3,1,3)
    plot(ca_15_resp)
end

%% lag curves at each voxel 2D (mean epoch)

if strcmp(mode_paradigm, 'rs')
    h = figure; set(h, 'WindowStyle', 'docked', 'NumberTitle', 'off', 'Name', ['lag 2D f1 meanEpoch'])
    for is=1:n_slices
        subplot(1,3,is)
        imagesc(lag*TR, 1:size(xc_f1_mean{is}, 2), xc_f1_mean{is}')
        cmin = min(min(xc_f1_mean{is}'));
        cmax = max(max(xc_f1_mean{is}'));
        title (['slice',num2str(is)]);
        xlabel('Times(s)');
        ylabel ('Voxels');
    end
    for is=1:n_slices
        subplot(1,3,is)
        colorbar
        caxis([cmin cmax])
    end
    
    h = figure; set(h, 'WindowStyle', 'docked', 'NumberTitle', 'off', 'Name', ['lag 2D f2 meanEpoch'])
    for is=1:n_slices
        subplot(1,3,is)
        imagesc(lag*TR, 1:size(xc_f2_mean{is}, 2), xc_f2_mean{is}')
        cmin = min(min(xc_f2_mean{is}'));
        cmax = max(max(xc_f2_mean{is}'));
        title (['slice',num2str(is)]);
        xlabel('Times(s)');
        ylabel ('Voxels');
    end
    for is=1:n_slices
        subplot(1,3,is)
        colorbar
        caxis([cmin cmax])
    end
    
    h = figure; set(h, 'WindowStyle', 'docked', 'NumberTitle', 'off', 'Name', ['lag 2D 15 meanEpoch'])
    for is=1:n_slices
        subplot(1,3,is)
        imagesc(lag*TR, 1:size(xc_15_mean{is}, 2), xc_15_mean{is}')
        cmin = min(min(xc_15_mean{is}'));
        cmax = max(max(xc_15_mean{is}'));
        title (['slice',num2str(is)]);
        xlabel('Times(s)');
        ylabel ('Voxels');
    end
    for is=1:n_slices
        subplot(1,3,is)
        colorbar
        caxis([cmin cmax])
    end
end

%% freq. band & mean fMRI slice plots
for ir = 1:n_plots
    % freq 1
    h = figure; set(h, 'WindowStyle', 'docked', 'NumberTitle', 'off', ...
        'Name', ['fmri ca-freq-15 s', num2str(is), ' r', num2str(runs(ir))])
    for is = 1:n_slices
        s_mean = mean(ctx{is}{ir}, 2);
        subplot(3,1,is)
        plot(t_fmri, zscore(s_mean)), hold on
        plot(t_fmri, zscore(ca_freq_15{ir}), 'LineWidth', 2)
        [cc, l] =xcorr(zscore(s_mean),zscore(ca_freq_15{ir}), 100, 'coeff');
        [~,ml] = max(cc);
        axis tight, box off
        title(['slice ', num2str(is),'   cc: ', num2str(cc(l==0)), '   max lag: ', num2str(l(ml))])
    end
    % freq 2
    %     h = figure; set(h, 'WindowStyle', 'docked', 'NumberTitle', 'off', ...
    %             'Name', ['fmri ca-freq-05 s', num2str(is), ' r', num2str(runs(ir))])
    %     for is = 1:n_slices
    %         s_mean = mean(ctx{is}{ir}, 2);
    %         subplot(3,1,is)
    %         plot(t_fmri, zscore(s_mean)), hold on
    %         plot(t_fmri, zscore(ca_freq_05{ir}), 'LineWidth', 2)
    %         [cc, l] =xcorr(zscore(s_mean),zscore(ca_freq_05{ir}), 100, 'coeff');
    %         axis tight, box off
    %         title(['slice ', num2str(is),'   cc: ', num2str(cc(l==0)), '   max lag: ', num2str(l(ml))])
    %     end
    % freq 3
    %     h = figure; set(h, 'WindowStyle', 'docked', 'NumberTitle', 'off', ...
    %             'Name', ['fmri ca-freq-01 s', num2str(is), ' r', num2str(runs(ir))])
    %     for is = 1:n_slices
    %         s_mean = mean(ctx{is}{ir}, 2);
    %         subplot(3,1,is)
    %         plot(t_fmri, zscore(s_mean)), hold on
    %         plot(t_fmri, zscore(ca_freq_01{ir}), 'LineWidth', 2)
    %         [cc, l] =xcorr(zscore(s_mean),zscore(ca_freq_01{ir}), 100, 'coeff');
    %         axis tight, box off
    %         title(['slice ', num2str(is),'   cc: ', num2str(cc(l==0)), '   max lag: ', num2str(l(ml))])
    %     end
end

%% filtered Ca & mean fMRI slice plots
for ir = 1:n_plots
    % filt 1
    h = figure; set(h, 'WindowStyle', 'docked', 'NumberTitle', 'off', ...
        'Name', ['fmri ca-f-001-01 s', num2str(is), ' r', num2str(runs(ir))])
    for is = 1:n_slices
        s_mean = mean(ctx{is}{ir}, 2);
        subplot(3,1,is)
        plot(t_fmri, zscore(s_mean)), hold on
        plot(t_fmri, zscore(ca_f1_ds{ir}), 'LineWidth', 2)
        [cc, l] =xcorr(zscore(s_mean),zscore(ca_f1_ds{ir}), 100, 'coeff');
        axis tight, box off
        title(['slice ', num2str(is),'   cc: ', num2str(cc(l==0)), '   max lag: ', num2str(l(ml))])
    end
    
    % filt 2
    h = figure; set(h, 'WindowStyle', 'docked', 'NumberTitle', 'off', ...
        'Name', ['fmri ca-f-001-20 s', num2str(is), ' r', num2str(runs(ir))])
    for is = 1:n_slices
        s_mean = mean(ctx{is}{ir}, 2);
        subplot(3,1,is)
        plot(t_fmri, zscore(s_mean)), hold on
        plot(t_fmri, zscore(ca_f2_ds{ir}), 'LineWidth', 2)
        [cc, l] =xcorr(zscore(s_mean),zscore(ca_f2_ds{ir}), 100, 'coeff');
        axis tight, box off
        title(['slice ', num2str(is),'   cc: ', num2str(cc(l==0)), '   max lag: ', num2str(l(ml))])
    end
    
    % filt 3
    %     h = figure; set(h, 'WindowStyle', 'docked', 'NumberTitle', 'off', ...
    %             'Name', ['fmri ca-f-0-002 s', num2str(is), ' r', num2str(runs(ir))])
    %     for is = 1:n_slices
    %         s_mean = mean(ctx{is}{ir}, 2);
    %         subplot(3,1,is)
    %         plot(t_fmri, zscore(s_mean)), hold on
    %         plot(t_fmri, zscore(ca_f3_ds{ir}), 'LineWidth', 2)
    %         [cc, l] =xcorr(zscore(s_mean),zscore(ca_f3_ds{ir}), 100, 'coeff');
    %         axis tight, box off
    %         title(['slice ', num2str(is),'   cc: ', num2str(cc(l==0)), '   max lag: ', num2str(l(ml))])
    %     end
end