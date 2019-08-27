params.complex = 1;
params.filter=[0,0,0,0];
params.phase = 0; 
params.zf = [0,0,0,0];
params.ft = [1,1,0,0];%2D 
font_size = 18;
axs_line_width = 2;
plot_linewidth = 2;
TR = 0.1; % one TR for one epoch
Fs = 1/TR; % Fs for sampling frequency
Spatial_Resolution = 0.05;
nCortical_Voxel = 45; 
offset_max = 2;

addpath('/Users/wangqi/Downloads/Matlab_scripts/Matlab_scripts/Bruker_Raw_Read')
addpath('/Users/wangqi/Downloads/Matlab_scripts/Matlab_scripts/utils/')
path = '/Users/wangqi/Downloads/From_Qi/12112018/37';
%% plug in data
path_fid{1,1}= strcat(path,'/','fid');
[image, protocol, pathstr] = BrReadImage_AutoRun(path_fid, params);
cd(pathstr);
nImage = size(image,3);
for num_image = 1 : nImage
                    close all;
                    if num_image == 1
                        cd(pathstr);
                        mkdir Results_Slice1;
                        pathstr1 = strcat(pathstr,'/Results_Slice1');
                        cd(pathstr1);
                    end
end
image_temp = squeeze(image(:,:,num_image,:,:)); % omit last 2 dimension
[nX, nY, nSlice] = size(squeeze(image_temp)); % nSlice is number of image inside per excitation, 0.1 sec per image
limit_WM = nX;
slice_num = 10;
% slice_num = 1;
%% display single linescanning line in k-sapce
imshow(abs(image_temp(:,:,slice_num)), []);% show single slice of raw linescanning image.
selected_kspace(:,:,:) = fft2c(image_temp(:,:,:));
h = figure();
imshow(abs(selected_kspace(:,:,slice_num)),[]);
title('fftline');% fft applied on single slice line scanning image.
set(gcf,'color','w');
set(gca,'FontSize',font_size,'FontWeight','Bold');
set(gca,'linewidth',axs_line_width);

imagesc(angle(selected_kspace(:,:,slice_num)));% imaged color map on angle of image after fft2c 
temp_line_profile(:,:,:) = ifft2c(selected_kspace(:,:,:));
temp0_cortical_depth_map(:,:,:) = fftc(temp_line_profile(:,:,:),2);% apply fftc on kspace images
colormap(jet);% down-sampling with more contrasted color
colorbar;
strName = sprintf('#%d Repeat Kspace(PhaseMap)',slice_num);
title(char(strName));
%% drawing a spectrum of 0-640 sec signal color map
for im_y = 1 : nY % 32 lines
    for slice_num = 1 : nSlice % 200 repeat
    temp_cortical_depth_map(:,im_y+nY*(slice_num-1)) = temp0_cortical_depth_map(:,im_y,slice_num);% 'im_y+nY*(Slice_num-1)'?
    end
end
total_temp_cortical_depth_map= temp_cortical_depth_map;
clear temp_cortical_depth_map;
h=figure();
imagesc(abs(total_temp_cortical_depth_map(:,:))); % from 0-640 sec, all response signal in sprectrum
% colormap(gray);
colormap(jet);
colorbar;
strName1 = sprintf('Cortical Depth');
title(strName1);
xlim([0 640*(1/TR)]) % set limit for x-axis from 0-640sec
xticks(0:160*(1/TR):640*(1/TR)) % set scale of x-axis from 0-640 sec
xticklabels({'0 sec','160 sec','320 sec','480s','640 sec'})
array_index4 = 0:0.5/Spatial_Resolution:size(total_temp_cortical_depth_map,1);
array_index4(1) = 1;
yticks(array_index4);
yticklabels({'0','0.5','1.0','1.5','2.0','2.5','3.0','3.5','4.0','4.5','5.0','5.5','6.0', ...
   '6.5','7.0','7.5','8.0','8.5','9.0','9.5','10.0','10.5','11.0','11.5','12.0','12.5' })
ylabel('Cortical depth (mm)');

slice_num=1;% Scan Number for Line Scanning

%% plot cortical surface.png
h=figure();
plot(abs(temp_line_profile(1:nX,nY/2+1,slice_num)),1:nX,'-'); 
set(gcf,'color','w');
set(gca,'FontSize',font_size,'FontWeight','Bold');
set(gca,'linewidth',axs_line_width);
strName2 = sprintf('Cortical Surface');
title(strName2);
set(gca, 'XAxisLocation', 'top');
set(gca,'Ydir','reverse');
%% Half max index location map
 for slice_num = 1 : size(total_temp_cortical_depth_map,2)
                        %         half_max_index = max(find(abs(temp_line_profile(1:max_i-1,nY/2+1,slice_num))< v/2))
                        %         [v, max_i]=max(abs(total_temp_cortical_depth_map(:,slice_num)),[],1);
    [v, max_i]=max(abs(total_temp_cortical_depth_map(1:limit_WM,slice_num)),[],1);
    if (isempty(find(abs(total_temp_cortical_depth_map(1:max_i-1,slice_num))< v/2, 1, 'last' )) || (find(abs(total_temp_cortical_depth_map(1:max_i-1,slice_num))< v/2, 1, 'last'))> (size(temp_line_profile,1)-nCortical_Voxel))
        if (slice_num >1 && temp_half_max_index(slice_num-1,1) ~= (size(temp_line_profile,1)-(nCortical_Voxel-1))) % added 02182019
            temp_half_max_index(slice_num,1) = temp_half_max_index(slice_num-1,1);
        else
            temp_half_max_index(slice_num,1) = size(temp_line_profile,1)-(nCortical_Voxel-1);
        end
    else
        temp_half_max_index(slice_num,1) = find(abs(total_temp_cortical_depth_map(1:max_i-1,slice_num))< v/2, 1, 'last' );
    end
end
figure, plot(abs(temp_half_max_index),'*b');
title('half max index location');
half_max_index = ceil(mean(temp_half_max_index,1))+ offset_max;
%% Cut off cortical map
total_cortical_depth_map = abs(total_temp_cortical_depth_map);
if (half_max_index+(nCortical_Voxel-1)) > size(total_cortical_depth_map,1)
    cortical_depth_map = total_cortical_depth_map(end-(nCortical_Voxel-1):end,:); % 50um spatial resolution, ~2mm
else
    cortical_depth_map = total_cortical_depth_map(half_max_index:half_max_index+(nCortical_Voxel-1),:); % 50um or 100um spatial resolution, ~2mm
end

[nX2, nY2] = size(cortical_depth_map);

h=figure();
imagesc(abs(cortical_depth_map(:,:)));
colormap(jet);
colorbar;
strName0 = 'Cutoff Cortical Depth';
title(strName0);
set(gcf,'color','w');
set(gca,'FontSize',font_size,'FontWeight','Bold');
set(gca,'linewidth',axs_line_width);
xlim([0 640*(1/TR)])
xticks(0:160*(1/TR):640*(1/TR))
xticklabels({'0s','160s','320s','480s','640s'})
%                 xlabel('Time(sec)');
array_index3 = 0:ceil(0.5/Spatial_Resolution):nCortical_Voxel;
array_index3(1) = 1;
yticks(array_index3)
yticklabels({'0','0.5','1.0','1.5','2.0','2.5','3.0'});
ylabel('Cortical depth (mm)')
%% Normalized cortical map[CHECK THIS BLOCK]
norm_cortical_map_temp = abs(cortical_depth_map(:,:))./max(max(abs(cortical_depth_map(:,:))));% ratio of each pixel to max pixel
Threshold_SI_MAP = abs(cortical_depth_map(:,:))./mean(norm_cortical_map_temp(:,:),2);% ratio of each pixel to average of rows 
norm_cortical_map = abs(Threshold_SI_MAP(:,:))./max(max(abs(Threshold_SI_MAP(:,:))));% idividual ratio to maxmium ratio of pixel

pre_stim = 10;
%nRepeat = 1;
Threshold_SI_percentage = zeros(size(cortical_depth_map,1),nSlice);% This var is single excitation
for nRepeat = 1 : nY

    SI_base_mean = mean(abs(cortical_depth_map(:,1+nSlice*(nRepeat-1):10+nSlice*(nRepeat-1))),2); % 1s off
    SI_i = abs(cortical_depth_map(:,1+nSlice*(nRepeat-1):nSlice*(nRepeat)));

    %     SI_percentage = (SI_i - SI_base_mean)./ SI_base_mean;
    for t = 1 : size(SI_i,2)
        SI_percentage(:,t) = (SI_i(:,t) - SI_base_mean(:,1))./ SI_base_mean(:,1);
    end
    %     SI_percentage(find((SI_percentage)<-0.01))=0; %hard thresolding
    %     SI_percentage(find(abs(SI_percentage)>0.1 | abs(SI_percentage)<0.01))=0; %hard thresolding

    Threshold_SI_percentage = Threshold_SI_percentage + SI_percentage;
end
Threshold_SI_percentage = Threshold_SI_percentage/nY;

data = Threshold_SI_percentage;
%// Define integer grid of coordinates for the above data
[X,Y] = meshgrid(1:size(data,2), 1:size(data,1));
%// Define a finer grid of points
[X2,Y2] = meshgrid(1:0.1:size(data,2), 1:0.1:size(data,1));
%// Interpolate the data and show the output
outData = interp2(X, Y, data, X2, Y2, 'Linear');
[nX2, nY2] = size(norm_cortical_map);
TA = nY2*TR;
h=figure();
max_val = 0.2;
imagesc(outData,[-0.01 max_val]);
imagesc((norm_cortical_map(:,:)));
%                 colormap(gray);
colormap(jet);
caxis([0.3 1])
colorbar;
title('Normalized Cortical Depth');
set(gcf,'color','w');
set(gca,'FontSize',font_size,'FontWeight','Bold');
set(gca,'linewidth',axs_line_width);
xlim([0 TA*(1/TR)])
xticks(0:TA/4*(1/TR):TA*(1/TR))
xticklabels({'0 sec','160 sec','320 sec','480 sec','640 sec'})
Ylim = get(gca, 'ylim');
set(gca, 'YTick', linspace(Ylim(1), Ylim(end), 10));
yticklabels({'0','','0.5','','1.0','','1.5','','2.0','2.25'})
ylabel('Cortical depth (mm)')
set(gcf, 'Position', [500, 500, 980, 300])
%% Averaged voxel with whole time course, using zscore method
plot(1:size(cortical_depth_map,2),zscore(mean(abs(cortical_depth_map(:,:)),1)),'-k','LineWidth',plot_linewidth);% lpot normalized signal w/ averaged columns
strName1 = sprintf('Averaged voxel with whole time course');
xlim([0 640*(1/TR)]);
xticks(0:160*(1/TR):640*(1/TR))
xticklabels({'0 sec','160 sec','320 sec','480 sec','640 sec'});
ylabel('fMRI Signal (a.u)');
ylim([-4 4])
title(strName1);
set(gcf,'color','w');
set(gcf, 'Position', [500, 500, 900, 300])
set(gca,'FontSize',font_size,'FontWeight','Bold');
set(gca,'linewidth',axs_line_width)
%% Averaged PSD with whole time course, using pwelch method
x = zscore(sum(abs(cortical_depth_map(:,:)),1)); % show trend of signal via z-score weighting
[pxx,f] = pwelch(x,[],[],[],Fs, 'onesided'); % transfer value of zscore to freq. domain, and do pwelch
h= figure();
plot(f,10*log10(pxx),'k','LineWidth',plot_linewidth);
ylabel('PSD');
ylim([-60 60])
xticklabels('manual')
xlim([0 5])
xticks(0:1:5)
xticklabels({'0','1','2','3','4','5'})
xlabel('Frequency (Hz)');
set(gcf, 'Position', [500, 500, 900, 300])
set(gca,'FontSize',font_size,'FontWeight','Bold');
set(gca,'linewidth',axs_line_width);
set(gcf,'color','w');
strName2 ='Averaged PSD with whole time course';
title(strName2);
%% Zoomed PSD w/ whole time course
norm_pxx = 10*log10(pxx)./ max(10*log10(pxx)); %what's meaning?
h= figure();
plot(f,norm_pxx,'k','LineWidth',plot_linewidth)
ylabel('Normalized PSD (a.u)');
xticklabels('manual')
xlim([0.01 0.1])
array_index8 = 0.01:0.01:0.1;
xticks(array_index8)
xticklabels({'0.01','0.02','0.03','0.04','0.05', ...
'0.06','0.07','0.08','0.09','0.1'})
 xlabel('Frequency (Hz)');
set(gcf, 'Position', [500, 500, 900, 300])
set(gca,'FontSize',font_size,'FontWeight','Bold');
set(gca,'linewidth',axs_line_width);
set(gcf,'color','w');
strName2 ='Zoomed PSD (0p01-0p1) with whole time course';
title(strName2);
%% Average time course
numberOfDataSets = nCortical_Voxel;
myColorOrder = jet(numberOfDataSets);
for temp = 1 : numberOfDataSets 
    plot(1:size(Threshold_SI_percentage,2), Threshold_SI_percentage(temp,:),'-');
    set(gca, 'ColorOrder', myColorOrder, 'NextPlot', 'replacechildren');
    hold on;
end
hold off
set(gcf,'color','w');
set(gca,'FontSize',font_size,'FontWeight','Bold');
set(gca,'linewidth',axs_line_width);
colormap(myColorOrder)
c =colorbar('YTickLabel',{'0', '2.25'},'YTick', [0 1], 'Ydir','reverse');
c.Label.String = 'Cortical depth (mm)';
c.Label.Rotation = 90;
set(c,'YTickMode','manual');
ylabel('Percentage (%)');
max_value = 0.3;
min_value = -0.05;

ylim([min_value max_value])
yticklabels('manual');
yticks(min_value:0.05:max_value)
yticklabels({'-5','0','5','10','15','20','25','30'});
xticklabels('manual');
xlim([0 200])
xticks(0:50:200)
xticklabels({'0 sec','5 sec','10 sec','15 sec','20 sec'})
strName1 = sprintf('Averaged Time Course');
title(strName1);

%% voxel-to-voxel correlation
norm_cortical_map = zscore(squeeze(abs(cortical_depth_map(:,:))),0,2)';

R = corrcoef(norm_cortical_map);% return correlation of 'norm_cortical_map'

Auto_Corr_Matrix_Setzero = R-diag(diag(R));% diag(diag) plot only diagnal element w/ rest on 0
corr_min = -0.2;
corr_max = 1.0;
h=figure();

Average_voxel_layers(1,:)  = mean(cortical_depth_map(1:3,:),1);%L1,0~0.15 mm, averaged voxel
Average_voxel_layers(2,:)  = mean(cortical_depth_map(4:11,:),1);%L2/3,0.15~0.55 mm, averaged voxel
Average_voxel_layers(3,:)  = mean(cortical_depth_map(12:16,:),1);%L4,0.55~0.8 mm, averaged voxel
Average_voxel_layers(4,:)  = mean(cortical_depth_map(17:27,:),1);%L5,0.8~1.4 mm, averaged voxel
Average_voxel_layers(5,:)  = mean(cortical_depth_map(28:40,:),1);%L6,1.4~2.0 mm, averaged voxel
Average_voxel_layers(6,:)  = mean(cortical_depth_map(41:nCortical_Voxel,:),1);%WM,2.1~3.0 mm, averaged voxel
imagesc(Auto_Corr_Matrix_Setzero);
colormap(jet);
c = colorbar('Location','eastoutside','Ticks',[corr_min corr_max]);
c.Label.String = 'Correlation Coefficient';
c.Label.Rotation = 270; % to rotate the text
strName1 = sprintf('Voxel-by-voxel Corrleation');

set(gca, 'XAxisLocation', 'top');
set(gcf,'color','w');
array_index3 = 0:ceil(0.25/Spatial_Resolution):nCortical_Voxel;
array_index3(1) = 1;
yticks(array_index3)
yticklabels({'0','','0.5','','1.0','','1.5','','2.0','2.25'})
ylabel('Cortical depth (mm)')
xticks(array_index3)
xticklabels({'0','','0.5','','1.0','','1.5','','2.0','2.25'})
xlabel('Cortical depth (mm)');
caxis([corr_min corr_max]);
axis image;
set(gca,'FontSize',font_size,'FontWeight','Bold');
set(gca,'linewidth',axs_line_width);