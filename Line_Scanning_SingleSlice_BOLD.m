clear;
close all;
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% made by S.Choi @MPI, Tuebingen
% initial version on 0501 2018
% modified : finding cortical surface index on 0912 2 018
% modified : bug fixed, miss matching fft and ifft with 'BrReadImage function' on 0927 2018
% modified : For PV6, protocol,NR == 32 added
% modified : multiple trials and animal with one click run on 0808 2019
params.complex = 1;
params.filter=[0,0,0,0];
params.phase = 0; 
params.zf = [0,0,0,0];
params.ft = [1,1,0,0];%2D 


addpath('/Users/levi/Downloads/Single_Slice_LINE_data_for_paper/Matlab_scripts/utils'); 
addpath('/Users/levi/Downloads/Single_Slice_LINE_data_for_paper/Matlab_scripts/Bruker_Raw_Read');  
% addpath('/Users/Daniel/2_Bruker_Scanner/Calcium_by_Filip');


% for biopac data
% physio_path = '/Users/Daniel/2_Bruker_Scanner/1_Bruker_Data/2018yrs/03122018/Biopac_Matfile/';
RegressPhysio_On = 0;

font_size = 18;
axs_line_width = 2;
plot_linewidth = 2;

% TR = 0.1;
TR = 0.05;
TR50  = 1;

% TR = 0.1;
Fs =1/TR;
% Spatial_Resolution = 0.1; %0.1
Spatial_Resolution = 0.05; %0.05
if Spatial_Resolution == 0.1
     nCortical_Voxel = 20;
elseif Spatial_Resolution == 0.05
%     nCortical_Voxel = 40; 
    nCortical_Voxel = 45; 
end
%%%%%%%%%%%%%WM
% limit_WM = 128;

offset_max = 2; % set conservative value
% offset_max = 4; %02212019
% Data Path
% path = '/Users/Daniel/2_Bruker_Scanner/1_Bruker_Data/2018yrs/03122018/';
for Animal = 4 : 5
    for RS = 0 : 0  % 1 : resting-state, 0 : task-related
        if RS == 0 %Task
             % Data Path
            switch Animal
                case 1
                    disp('Animal #1')
                    %set your data path in next line
                    path = '/Users/levi/Downloads/Single_Slice_LINE_data_for_paper/03122018/';
                    Exp_Number_Array1 = [63 65 67 69]; %03122018
                case 2
                    disp('Animal #2')
                     %set your data path in next line
                    path = '/Users/levi/Downloads/Single_Slice_LINE_data_for_paper/05062018/'; 
                    Exp_Number_Array1 = [38 40 42]; %05062018
                case 3
                    disp('Animal #3')
                     %set your data path in next line
                    path = '/Users/levi/Downloads/Single_Slice_LINE_data_for_paper/07052018/'; 
                    Exp_Number_Array1 = [25 27 37 39]; %07052018
                case 4
                    disp('Animal #4')
                    %set your data path in next line
                    path = '/Users/levi/Downloads/Single_Slice_LINE_data_for_paper/08162018/'; 
                    if TR50 == 1
                        Exp_Number_Array1 = [37 40 41 42]; %08162018
                    else
                        Exp_Number_Array1 = [31 32 36]; %08162018
                    end
                case 5
                    disp('Animal #5')
                     %set your data path in next line
                    path = '/Users/levi/Downloads/Single_Slice_LINE_data_for_paper/09252018/'; 
                    if TR50 == 1
                        Exp_Number_Array1 = [29 31 34]; %08162018
                    else
                        Exp_Number_Array1 = [18 23 25 26 28 ]; %09252018
                    end
                otherwise
                    disp('NO Input!')
                    exit
            end
        else %% Resting
            % Data Path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% take out path of each trail w.r.t. animals %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            switch Animal
                case 1
                    disp('Animal #1')
                    %set your data path in next line
                    path = '/Users/Daniel/2_Bruker_Scanner/1_Bruker_Data/2018yrs/03122018/'; 
                    Exp_Number_Array1 = [64 66 68 70]; %03122018
                case 2
                    disp('Animal #2')
                     %set your data path in next line
                    path = '/Users/Daniel/2_Bruker_Scanner/1_Bruker_Data/2018yrs/05062018/'; 
                    Exp_Number_Array1 = [39 41]; %05062018
                case 3
                    disp('Animal #3')
                     %set your data path in next line
                    path = '/Users/Daniel/2_Bruker_Scanner/1_Bruker_Data/2018yrs/07052018/'; 
                    Exp_Number_Array1 = [28 38]; %07052018
                case 4
                    disp('Animal #4')
                    %set your data path in next line
                    path = '/Users/Daniel/2_Bruker_Scanner/1_Bruker_Data/2018yrs/08162018/'; 
                    Exp_Number_Array1 = [33 34 39]; %08162018
                case 5
                    disp('Animal #5')
                     %set your data path in next line
                    path = '/Users/Daniel/2_Bruker_Scanner/1_Bruker_Data/2018yrs/09252018/'; 
                    Exp_Number_Array1 = [22 24 27 30 32 33 35]; %09252018
                otherwise
                    disp('NO Input!')
                    exit
            end
           
            %          Exp_Number_Array1 = [33]; %02072019
            
        end
        
        for exp_case_num = 1 : size(Exp_Number_Array1,2)  %read how many columns in 'Exp_Number_Array1', which is array of trails
            close all;
            
            data_num2str = num2str(Exp_Number_Array1(exp_case_num));
            path_temp = strcat('Test #',data_num2str)
            path_fid{1,1}= strcat(path,'/',data_num2str,'/fid');
            [image, protocol, pathstr] = BrReadImage_AutoRun(path_fid, params);
            % strPath = 30;
            % newPath = strcat(path, num2str(strPath));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create new folder for results of each slices, w/wo filter %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            cd(pathstr);
            nImage = size(image,3);
            %          nImage = 2;
            for Filter_On = 0 : 1
                close all
                for num_image = 1 : nImage
                    %             for num_image = 2 : 2
                    close all;
                    if num_image == 1
                        cd(pathstr);
                        mkdir Results_Slice1;
                        pathstr1 = strcat(pathstr,'/Results_Slice1');
                        cd(pathstr1);
                    elseif num_image == 2
                        close all;
                        cd(pathstr);
                        mkdir Results_Slice2;
                        pathstr1 = strcat(pathstr,'/Results_Slice2');
                        cd(pathstr1);
                    elseif num_image == 3
                        close all;
                        cd(pathstr);
                        mkdir Results_Slice3;
                        pathstr1 = strcat(pathstr,'/Results_Slice3');
                        cd(pathstr1);
                    end
                    image_temp = squeeze(image(:,:,num_image,:,:)); % Unsolved sentence 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reset bandpass filter with 0.01-0.1, and filtering %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if Filter_On == 1
                        mkdir Results_Filtering_0p1;
                        pathstr2= strcat(pathstr1,'/Results_Filtering_0p1');
                        cd(pathstr2);
                        Fs = 1/TR;
                        %                     order    = 2048; %use it now
                        order    = 4096; %use it now
                        
                        fcutlow  = 0.01;   %low cut frequency in Hz 0.01 %play with this number matters a lot
                        fcuthigh = 0.1;   %intend 0.1 Hz but consider transient band
                        Fc1 = fcutlow;
                        Fc2 = fcuthigh;
                        
                        N = order;
                        beta = 0.005;
                        win = kaiser(N+1, beta);
                        flag = 'scale';  % Sampling Flag
                        
                        % Calculate the coefficients using the FIR1 function.
                        b  = fir1(N, [Fc1, Fc2]/(Fs/2), 'bandpass', win, flag);
                        
                        fvtool(b,1,'Fs',Fs)
                        Hd = dfilt.dffir(b);
                    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forget about RegressPhysio_On first  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if RegressPhysio_On == 1
                        if Filter_On == 0
                            mkdir Results_RegressPhysio;
                            pathstr3 = strcat(pathstr1,'/Results_RegressPhysio');
                            cd(pathstr3);
                        elseif Filter_On == 1
                            mkdir Results_Filter0p1_RegressPhysio;
                            pathstr2= strcat(pathstr1,'/Results_Filtering_0p1');
                            pathstr3 = strcat(pathstr2,'/Results_Filter0p1_RegressPhysio');
                            cd(pathstr3);
                        end
                    end
                    [nX, nY, nSlice] = size(squeeze(image_temp));
                    limit_WM = nX;
                    if Animal == 1
                        nSlice =200;
                    end
                    squeezed_image = squeeze(image_temp); % delete the dimension where value is 1
                    clear selected_kspace;
                    slice_num = 10; %center slice
                    figure();
                    imshow(abs(squeezed_image(:,:,slice_num)), []);
                    selected_kspace(:,:,:) = fft2c(squeezed_image(:,:,:));
                    
                    h=figure();
                    imshow((abs(selected_kspace(:,:,slice_num))), []);
                    set(gcf,'color','w');
                    set(gca,'FontSize',font_size,'FontWeight','Bold');
                    set(gca,'linewidth',axs_line_width);
                    strName = sprintf('#%d Repeat Kspace',slice_num);
                    title(char(strName));
                    saveas(h,strName,'png');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generating K-space images, which is a line%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                  
                    h=figure();
                    imagesc(((abs(selected_kspace(:,:,slice_num)))));
                    set(gcf,'color','w');
                    set(gca,'FontSize',font_size,'FontWeight','Bold');
                    set(gca,'linewidth',axs_line_width);
                    strName = sprintf('#%d Repeat Kspace(Log-scale)',slice_num);
                    title(char(strName));
                    saveas(h,strName,'png');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generating K-space images, which is a line%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    h=figure();
                    imagesc(((angle(selected_kspace(:,:,slice_num)))));
                    colormap(jet);
                    colorbar;
                    set(gcf,'color','w');
                    set(gca,'FontSize',font_size,'FontWeight','Bold');
                    set(gca,'linewidth',axs_line_width);
                    strName = sprintf('#%d Repeat Kspace(PhaseMap)',slice_num);
                    title(char(strName));
                    saveas(h,strName,'png');
                    
                    clear temp_line_profile;
                    temp_line_profile(:,:,:) = ifft2c(selected_kspace(:,:,:)); % 2D fft
                    h=figure();
                    imshow(abs(temp_line_profile(:,:,slice_num)), []);
                    strName1 = sprintf('#%d Line Profile',slice_num);
                    title(strName1);
                    set(gcf,'color','w');
                    set(gca,'FontSize',font_size,'FontWeight','Bold');
                    set(gca,'linewidth',axs_line_width);
                    saveas(h,strName1,'png');
                    clear selected_kspace;
                    %%
                    clear temp0_cortical_depth_map;
                    temp0_cortical_depth_map(:,:,:) = fftc(temp_line_profile(:,:,:),2); % ifft along y direction, image to freq
                    clear temp_cortical_depth_map;
                    
                    for im_y = 1 : nY % 32 lines
                        for slice_num = 1 : nSlice % 200 repeat
                            temp_cortical_depth_map(:,im_y+nY*(slice_num-1)) = temp0_cortical_depth_map(:,im_y,slice_num);
                        end
                    end
                    total_temp_cortical_depth_map= temp_cortical_depth_map;
                    clear temp_cortical_depth_map;
                    h=figure();
                    imagesc(abs(total_temp_cortical_depth_map(:,:)));
                    colormap(gray);
                    strName1 = sprintf('Cortical Depth');
                    title(strName1);
                    set(gcf,'color','w');
                    set(gca,'FontSize',font_size,'FontWeight','Bold');
                    set(gca,'linewidth',axs_line_width);
                    xlim([0 640*(1/TR)])
                    xticks(0:160*(1/TR):640*(1/TR))
                    xticklabels({'0 sec','160 sec','320 sec','480s','640 sec'})
                    %                 xlabel('Time(sec)');
                    %     array_index4 = 0:10:120;
                    array_index4 = 0:0.5/Spatial_Resolution:size(total_temp_cortical_depth_map,1);
                    array_index4(1) = 1;
                    yticks(array_index4)
                    yticklabels({'0','0.5','1.0','1.5','2.0','2.5','3.0','3.5','4.0','4.5','5.0','5.5','6.0', ...
                        '6.5','7.0','7.5','8.0','8.5','9.0','9.5','10.0','10.5','11.0','11.5','12.0','12.5' })
                    ylabel('Cortical depth (mm)')
                    set(gcf, 'Position', [500, 500, 900, 600])
                    saveas(h,strName1,'png');
                    %%
                    slice_num=1;% Scan Number for Line Scanning
                    h=figure();
                    plot(abs(temp_line_profile(1:nX,nY/2+1,slice_num)),1:nX,'-');
                    set(gcf,'color','w');
                    set(gca,'FontSize',font_size,'FontWeight','Bold');
                    set(gca,'linewidth',axs_line_width);
                    strName2 = sprintf('Cortical Surface');
                    title(strName2);
                    set(gca, 'XAxisLocation', 'top');
                    set(gca,'Ydir','reverse');
                    saveas(h,strName2,'png');
                    
                    
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
                    
                    half_max_index = ceil(mean(temp_half_max_index,1))+ offset_max
                    
                    total_cortical_depth_map = abs(total_temp_cortical_depth_map);
                    clear total_temp_cortical_depth_map
                    
                    
                    
                    %%
                    clear cortical_depth_map;
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
                    %                 set(gcf, 'Position', [500, 500, 900, 600])
                    saveas(h,strName0,'png');
                    save cortical_depth_map cortical_depth_map;
                    
                    %%
                    if Filter_On > 0
                        clear temp_filtered_cortical_depth_map1
                        clear temp_filtered_cortical_depth_map2
                        %                     clear filtered_cortical_depth_map_new
                        origin_cortical_depth_map  = (cortical_depth_map);
                        filtered_cortical_depth_map = zeros(size(cortical_depth_map));
                        nVoxel = size(filtered_cortical_depth_map,1);
                        nTime = size(filtered_cortical_depth_map,2);
                        for voxel_num = 1 : nVoxel
                            
                            %prevent the signal fluctuation in front and behind part after filtering
                            temp_filtered_cortical_depth_map1(voxel_num,:)=fliplr((cortical_depth_map(voxel_num,2:N+1)));  % maintain continuity in level and slope
                            temp_filtered_cortical_depth_map2(voxel_num,:)=fliplr((cortical_depth_map(voxel_num,end-N:end-1)));
                            
                            filtered_cortical_depth_map_new(voxel_num,1:N)                 = temp_filtered_cortical_depth_map1(voxel_num,1:N);
                            filtered_cortical_depth_map_new(voxel_num,N+1:nTime+N)         = (cortical_depth_map(voxel_num,1:nTime));
                            filtered_cortical_depth_map_new(voxel_num,nTime+N+1:nTime+2*N) = temp_filtered_cortical_depth_map2(voxel_num,1:N);
                            
                            %                         filtered_cortical_depth_map(voxel_num,:) =  filter(Hd,squeeze(cortical_depth_map(voxel_num,:)));
                            filtered_cortical_depth_map_new(voxel_num,:) =  filter(Hd,squeeze(filtered_cortical_depth_map_new(voxel_num,:)));
                            delay = mean(grpdelay(Hd));
                            filtered_cortical_depth_map_new(voxel_num,:) = circshift(filtered_cortical_depth_map_new(voxel_num,:),(-1)*delay,2);
                        end
                        
                        %Cut unnessary part
                        filtered_cortical_depth_map(:,1:nTime) =  filtered_cortical_depth_map_new(:,N+1:N+nTime);
                        
                    else %Filter_Off
                        filtered_cortical_depth_map = cortical_depth_map; %without filtering
                    end
                    
                    %replace
                    clear cortical_depth_map;
                    cortical_depth_map = abs(filtered_cortical_depth_map);
                    
                    %%
                    
                    if Filter_On == 0 || Filter_On == 1
                        if num_image == 1
                            save half_max_index_s1 half_max_index;
                            save line_scanning_data_s1 total_cortical_depth_map;
                            save cortical_depth_map_s1 cortical_depth_map;
                        elseif num_image == 2
                            save half_max_index_s2 half_max_index;
                            save line_scanning_data_s2 total_cortical_depth_map;
                            save cortical_depth_map_s2 cortical_depth_map;
                        elseif num_image == 3
                            save half_max_index_s3 half_max_index;
                            save line_scanning_data_s3 total_cortical_depth_map;
                            save cortical_depth_map_s3 cortical_depth_map;
                            offset_index0 = 15;
                            nX_point0 = half_max_index + offset_index0;%L4 layaer
                            synchorized_single_voxel = total_cortical_depth_map(nX_point0,:);
                            synchorized_averaged_voxel = mean(abs(total_cortical_depth_map(half_max_index+1:half_max_index+nCortical_Voxel,:)),1);
                            save synchorized_single_voxel synchorized_single_voxel
                            save synchorized_averaged_voxel synchorized_averaged_voxel
                        end
                    end
                    
                    %%
                    if RS == 1 %resting-state
                        norm_cortical_map = abs(cortical_depth_map(:,:))./max(max(abs(cortical_depth_map(:,:))));
                    else%Task
                        pre_stim = 10;
                        %                 nRepeat = 1;
                        for nRepeat = 1 : nY
                            
                            SI_base_mean = mean(abs(cortical_depth_map(:,1+nSlice*(nRepeat-1):pre_stim+nSlice*(nRepeat-1))),2); % 1s off
                            SI_i = abs(cortical_depth_map(:,1+nSlice*(nRepeat-1):nSlice*(nRepeat)));
                            
                            %     SI_percentage = (SI_i - SI_base_mean)./ SI_base_mean;
                            for t = 1 : size(SI_i,2)
                                SI_percentage(:,t) = (SI_i(:,t) - (SI_base_mean(:,1)))./(SI_base_mean(:,1));
                            end
                            %     SI_percentage(find((SI_percentage)<-0.01))=0; %hard thresolding
                            %     SI_percentage(find(abs(SI_percentage)>0.1 | abs(SI_percentage)<0.01))=0; %hard thresolding
                            
                            Threshold_SI_MAP(:,1+nSlice*(nRepeat-1):nSlice*(nRepeat)) = SI_percentage;
                        end
                        %                 SI_base_mean = mean(abs(cortical_depth_map(:,1+nSlice*(nRepeat-1):pre_stim+nSlice*(nRepeat-1))),2); % 1s off
                        
                        %                 norm_cortical_map = Threshold_SI_MAP;
                        norm_cortical_map = abs(Threshold_SI_MAP(:,:))./max(max(abs(Threshold_SI_MAP(:,:))));
                    end
                    %                  norm_cortical_map = (Threshold_SI_MAP(:,:));
                    [nX2, nY2] = size(norm_cortical_map);
                    TA = nY2*TR;
                    h=figure();
                    imagesc((norm_cortical_map(:,:)));
                    %                 colormap(gray);
                    colormap(jet);
                    
                    if RS == 1
                        caxis([0.5 1])
                        colorbar;
                    else
                        
                        if Filter_On == 0
                            max_val = 0.1;
                            min_val = -0.1;
                            caxis([min_val max_val])
                            colorbar('YTickLabel',{strcat(num2str(min_val*100),'%'),strcat(num2str(max_val*100),'%')},'YTick', [min_val max_val]);
                        else
                            max_val = 0.3;
                            min_val = 0;
                            caxis([min_val max_val])
                            colorbar('YTickLabel',{strcat(num2str(min_val*100),'%'),strcat(num2str(max_val*100),'%')},'YTick', [min_val max_val]);
                        end
                    end
                    if num_image == 1
                        strName0 = sprintf('Percentage Cortical Depth');
                    else
                        strName0 = sprintf('Slice #%d, Percentage Cortical Depth',num_image);
                    end
                    title(strName0);
                    set(gcf,'color','w');
                    set(gca,'FontSize',font_size,'FontWeight','Bold');
                    set(gca,'linewidth',axs_line_width);
                    
                    xlim([0 TA*(1/TR)])
                    if TA == 640
                        xticks(0:TA/4*(1/TR):TA*(1/TR))
                        xticklabels({'0 sec','160 sec','320 sec','480 sec','640 sec'})
                    elseif TA == 768
                        xticks(0:192*(1/TR):TA*(1/TR))
                        xticklabels({'0 sec','192 sec','384 sec','576 sec','768 sec'})
                    else
                        disp('No match with this TA!')
                    end
                    if Spatial_Resolution == 0.05
                        if nCortical_Voxel == 60
                            Ylim = get(gca, 'ylim');
                            set(gca, 'YTick', linspace(Ylim(1), Ylim(end), 6));
                            yticklabels({'0','0.5','1.0','1.5','2.0','2.5','3.0'})
                            ylabel('Cortical depth (mm)')
                        elseif nCortical_Voxel == 45
                            Ylim = get(gca, 'ylim');
                            set(gca, 'YTick', linspace(Ylim(1), Ylim(end), 10));
                            yticklabels({'0','','0.5','','1.0','','1.5','','2.0','2.25'})
                            ylabel('Cortical depth (mm)')
                        else
                            Ylim = get(gca, 'ylim');
                            set(gca, 'YTick', linspace(Ylim(1), Ylim(end), 4));
                            yticklabels({'0','0.5','1.0','1.5','2.0','2.5','3.0'})
                            ylabel('Cortical depth (mm)')
                        end
                        
                    elseif Spatial_Resolution == 0.1
                        if nCortical_Voxel == 30
                            Ylim = get(gca, 'ylim');
                            set(gca, 'YTick', linspace(Ylim(1), Ylim(end), 6));
                        else
                            Ylim = get(gca, 'ylim');
                            set(gca, 'YTick', linspace(Ylim(1), Ylim(end), 4));
                        end
                        yticklabels({'0','0.5','1.0','1.5','2.0','2.5','3.0'})
                        ylabel('Cortical depth (mm)')
                    end
                    set(gcf, 'Position', [500, 500, 980, 300])
                    %             set(gcf, 'Position', [500, 500, 2500, 500])
                    saveas(h,strName0,'png');
                    
                    %%
                    norm_cortical_map_temp = abs(cortical_depth_map(:,:))./max(max(abs(cortical_depth_map(:,:))));
                    Threshold_SI_MAP = abs(cortical_depth_map(:,:))./mean(norm_cortical_map_temp(:,:),2);
                    norm_cortical_map = abs(Threshold_SI_MAP(:,:))./max(max(abs(Threshold_SI_MAP(:,:))));
                    %                  norm_cortical_map = (abs(Threshold_SI_MAP(:,:))-min(min(abs(Threshold_SI_MAP(:,:)))))./(max(max(abs(Threshold_SI_MAP(:,:))))--min(min(abs(Threshold_SI_MAP(:,:)))));
                    %      else%Task
                    %                  norm_cortical_map = (Threshold_SI_MAP(:,:));
                    [nX2, nY2] = size(norm_cortical_map);
                    TA = nY2*TR;
                    h=figure();
                    imagesc((norm_cortical_map(:,:)));
                    %                 colormap(gray);
                    colormap(jet);
                    
                    if RS == 1
                        caxis([0.5 1])
                        colorbar;
                    else
                        if TR50 == 1
                            caxis([0.4 1])
                        else
                            caxis([0.3 1])
                        end
                        %                     caxis([0.5 1])
                        %                     caxis([0.2 0.8])
                        colorbar;
                    end
                    %                 if num_image == 1
                    strName0 = sprintf('Normalized Cortical Depth');
                    %                 else
                    %                     strName0 = sprintf('Slice #%d, Normalized Cortical Depth',num_image);
                    %                 end
                    title(strName0);
                    set(gcf,'color','w');
                    set(gca,'FontSize',font_size,'FontWeight','Bold');
                    set(gca,'linewidth',axs_line_width);
                    
                    xlim([0 TA*(1/TR)])
                    if TA == 640
                        xticks(0:TA/4*(1/TR):TA*(1/TR))
                        xticklabels({'0 sec','160 sec','320 sec','480 sec','640 sec'})
                    elseif TA == 768
                        xticks(0:192*(1/TR):TA*(1/TR))
                        xticklabels({'0 sec','192 sec','384 sec','576 sec','768 sec'})
                    else
                        disp('No match with this TA!')
                    end
                    %                 array_index3 = 0:ceil(0.5/Spatial_Resolution):nCortical_Voxel;
                    %                 array_index3(1) = 1;
                    %                 yticks(array_index3)
                    %                 yticklabels({'0','0.5','1.0','1.5','2.0','2.5','3.0'});
                    %                 ylabel('Cortical depth (mm)')
                    if Spatial_Resolution == 0.05
                        if nCortical_Voxel == 60
                            Ylim = get(gca, 'ylim');
                            set(gca, 'YTick', linspace(Ylim(1), Ylim(end), 6));
                            yticklabels({'0','0.5','1.0','1.5','2.0','2.5','3.0'})
                            ylabel('Cortical depth (mm)')
                        elseif nCortical_Voxel == 45
                            Ylim = get(gca, 'ylim');
                            set(gca, 'YTick', linspace(Ylim(1), Ylim(end), 10));
                            yticklabels({'0','','0.5','','1.0','','1.5','','2.0','2.25'})
                            ylabel('Cortical depth (mm)')
                        else
                            Ylim = get(gca, 'ylim');
                            set(gca, 'YTick', linspace(Ylim(1), Ylim(end), 4));
                            yticklabels({'0','0.5','1.0','1.5','2.0','2.5','3.0'})
                            ylabel('Cortical depth (mm)')
                        end
                        
                    elseif Spatial_Resolution == 0.1
                        if nCortical_Voxel == 30
                            Ylim = get(gca, 'ylim');
                            set(gca, 'YTick', linspace(Ylim(1), Ylim(end), 6));
                        else
                            Ylim = get(gca, 'ylim');
                            set(gca, 'YTick', linspace(Ylim(1), Ylim(end), 4));
                        end
                        yticklabels({'0','0.5','1.0','1.5','2.0','2.5','3.0'})
                        ylabel('Cortical depth (mm)')
                    end
                    set(gcf, 'Position', [500, 500, 980, 300])
                    %             set(gcf, 'Position', [500, 500, 2500, 500])
                    saveas(h,strName0,'png');
                    
                    %%  diff_cortex = abs(cortex{1}{exp_case_num}')- abs(cortical_depth_map(:,:));
                    if RegressPhysio_On == 1
                        norm_cortical_map_diff = abs(diff_cortex)./max(max(abs(diff_cortex)));
                        %                  Threshold_SI_MAP = abs(cortex{1}{exp_case_num}')./mean(norm_cortical_map_temp(:,:),2);
                        %                  norm_cortical_map_orig = abs(Threshold_SI_MAP(:,:))./max(max(abs(Threshold_SI_MAP(:,:))));
                        
                        %                 diff_cortex = abs(abs(norm_cortical_map_orig)- abs(norm_cortical_map))./max(max(abs(abs(norm_cortical_map_orig)- abs(norm_cortical_map))));
                        
                        %                  norm_cortical_map = (abs(Threshold_SI_MAP(:,:))-min(min(abs(Threshold_SI_MAP(:,:)))))./(max(max(abs(Threshold_SI_MAP(:,:))))--min(min(abs(Threshold_SI_MAP(:,:)))));
                        %      else%Task
                        %                  norm_cortical_map = (Threshold_SI_MAP(:,:));
                        [nX2, nY2] = size(norm_cortical_map_diff);
                        TA = nY2*TR;
                        h=figure();
                        imagesc((norm_cortical_map_diff(:,:)));
                        %                 colormap(gray);
                        colormap(jet);
                        
                        if RS == 1
                            caxis([0.5 1])
                            colorbar;
                        else
                            caxis([0.1 1])
                            %                     caxis([0.5 1])
                            %                     caxis([0.2 0.8])
                            colorbar;
                        end
                        %                 if num_image == 1
                        strName0 = sprintf('Diff Normalized Cortical Depth');
                        %                 else
                        %                     strName0 = sprintf('Slice #%d, Normalized Cortical Depth',num_image);
                        %                 end
                        title(strName0);
                        set(gcf,'color','w');
                        set(gca,'FontSize',font_size,'FontWeight','Bold');
                        set(gca,'linewidth',axs_line_width);
                        
                        xlim([0 TA*(1/TR)])
                        if TA == 640
                            xticks(0:TA/4*(1/TR):TA*(1/TR))
                            xticklabels({'0 sec','160 sec','320 sec','480 sec','640 sec'})
                        elseif TA == 768
                            xticks(0:192*(1/TR):TA*(1/TR))
                            xticklabels({'0 sec','192 sec','384 sec','576 sec','768 sec'})
                        else
                            disp('No match with this TA!')
                        end
                        array_index3 = 0:ceil(0.5/Spatial_Resolution):nCortical_Voxel;
                        array_index3(1) = 1;
                        yticks(array_index3)
                        yticklabels({'0','0.5','1.0','1.5','2.0','2.5','3.0'});
                        ylabel('Cortical depth (mm)')
                        set(gcf, 'Position', [500, 500, 980, 300])
                        %             set(gcf, 'Position', [500, 500, 2500, 500])
                        saveas(h,strName0,'png');
                    end
                    %%
                    
                    if Filter_On > 0
                        h=figure();
                        imagesc(abs(cortical_depth_map(:,:)));
                        colormap(gray);
                        strName0 = sprintf('Filtered Cortical Depth');
                        title(strName0);
                        set(gcf,'color','w');
                        set(gca,'FontSize',font_size,'FontWeight','Bold');
                        set(gca,'linewidth',axs_line_width);
                        xlim([0 640*(1/TR)])
                        xticks(0:160*(1/TR):640*(1/TR))
                        xticklabels({'0 sec','160 sec','320 sec','480s','640 sec'})
                        %                     xlabel('Time(sec)');
                        %     array_index4 = 0:10:120;
                        array_index4 = 0:ceil(0.5/Spatial_Resolution):size(cortical_depth_map,1);
                        array_index4(1) = 1;
                        yticks(array_index4)
                        yticklabels({'0','0.5','1.0','1.5','2.0','2.5','3.0','3.5','4.0','4.5','5.0','5.5','6.0', ...
                            '6.5','7.0','7.5','8.0','8.5','9.0','9.5','10.0','10.5','11.0','11.5','12.0','12.5' })
                        ylabel('Filtered Cortical depth (mm)')
                        saveas(h,strName0,'png');
                    end
                    
                    %%
                    h=figure();
                    %     plot(1:size(synchorized_averaged_voxel,2),zscore(synchorized_averaged_voxel(1,:)),'-k','LineWidth',1.2);
                    %
                    plot(1:size(cortical_depth_map,2),zscore(mean(abs(cortical_depth_map(:,:)),1)),'-k','LineWidth',plot_linewidth);
                    %                 strName1 = sprintf('Slice #%d, Averaged voxel with whole time course',num_image);
                    strName1 = sprintf('Averaged voxel with whole time course');
                    xlim([0 640*(1/TR)])
                    xticks(0:160*(1/TR):640*(1/TR))
                    xticklabels({'0 sec','160 sec','320 sec','480 sec','640 sec'})
                    ylabel('fMRI Signal (a.u)');
                    ylim([-4 4])
                    title(strName1);
                    set(gcf,'color','w');
                    set(gcf, 'Position', [500, 500, 900, 300])
                    set(gca,'FontSize',font_size,'FontWeight','Bold');
                    set(gca,'linewidth',axs_line_width)
                    saveas(h,strName1,'png');
                    
                    
                    
                    %%
                    x = zscore(sum(abs(cortical_depth_map(:,:)),1));
                    [pxx,f] = pwelch(x,[],[],[],Fs, 'onesided');
                    h= figure();
                    plot(f,10*log10(pxx),'k','LineWidth',plot_linewidth)
                    ylabel('PSD');
                    ylim([-60 60])
                    xticklabels('manual')
                    if Fs == 10
                        xlim([0 5])
                        xticks(0:1:5)
                        xticklabels({'0','1','2','3','4','5'})
                    elseif Fs == 20
                        xlim([0 10])
                        xticks(0:1:10)
                        xticklabels({'0','1','2','3','4','5','6','7','8','9','10'})
                    end
                    xlabel('Frequency (Hz)');
                    set(gcf, 'Position', [500, 500, 900, 300])
                    set(gca,'FontSize',font_size,'FontWeight','Bold');
                    set(gca,'linewidth',axs_line_width);
                    set(gcf,'color','w');
                    strName2 ='Averaged PSD with whole time course';
                    title(strName2);
                    saveas(h,strName2,'png');
                    
                    %%
                    norm_pxx = 10*log10(pxx)./ max(10*log10(pxx));
                    %                                 norm_pxx = 10*log10(pxx);
                    h= figure();
                    plot(f,norm_pxx,'k','LineWidth',plot_linewidth)
                    ylabel('Normalized PSD (a.u)');
                    xticklabels('manual')
                    if Fs == 10
                        xlim([0.01 0.1])
                        array_index8 = 0.01:0.01:0.1;
                        xticks(array_index8)
                        xticklabels({'0.01','0.02','0.03','0.04','0.05', ...
                            '0.06','0.07','0.08','0.09','0.1'})
                    elseif Fs == 20
                        xlim([0.01 0.1])
                        xticks(0.01:0.01:0.1)
                        xticklabels({'0.01','0.02','0.03','0.04','0.05', ...
                            '0.06','0.07','0.08','0.09','0.1'})
                    end
                    xlabel('Frequency (Hz)');
                    set(gcf, 'Position', [500, 500, 900, 300])
                    set(gca,'FontSize',font_size,'FontWeight','Bold');
                    set(gca,'linewidth',axs_line_width);
                    set(gcf,'color','w');
                    strName2 ='Zoomed PSD (0p01-0p1) with whole time course';
                    title(strName2);
                    saveas(h,strName2,'png');
                    
                    
                    %%
                    
                    Threshold_SI_percentage = zeros(size(cortical_depth_map,1),nSlice);
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
                    Threshold_SI_percentage = Threshold_SI_percentage/nY; %for all time average not for single
                    
                    data = Threshold_SI_percentage;
                    %// Define integer grid of coordinates for the above data
                    [X,Y] = meshgrid(1:size(data,2), 1:size(data,1));
                    %// Define a finer grid of points
                    [X2,Y2] = meshgrid(1:0.1:size(data,2), 1:0.1:size(data,1));
                    %// Interpolate the data and show the output
                    outData = interp2(X, Y, data, X2, Y2, 'Linear');
                    
                    h=figure();
                    if Filter_On == 0
                        max_val = 0.2;
                        %                     max_val = 0.15;
                        %                     max_val = 0.1;
                        imagesc(outData,[-0.01 max_val]);
                    else
                        max_val = 0.5;
                        imagesc(outData,[-0.01 max_val]);
                    end
                    
                    strName105 = sprintf('Averaged Percentage Change');
                    
                    title(strName105);
                    colormap(jet);
                    if Filter_On == 0
                        colorbar('YTickLabel',{'-1%',strcat(num2str(max_val*100),'%')},'YTick', [-0.01 max_val]);
                    else
                        colorbar('YTickLabel',{'-1%',strcat(num2str(max_val*100),'%')},'YTick', [-0.01 max_val]);
                    end
                    
                    array_index6 = 0:size(outData,2)/4:size(outData,2);
                    array_index6(1) = 1;
                    xticks(array_index6)
                    xticklabels({'0 sec','5 sec','10 sec','15 sec','20 sec'})
                    %                 xlabel('Time(sec)');
                    if Spatial_Resolution == 0.05
                        if nCortical_Voxel == 60
                            Ylim = get(gca, 'ylim');
                            set(gca, 'YTick', linspace(Ylim(1), Ylim(end), 6));
                            yticklabels({'0','0.5','1.0','1.5','2.0','2.5','3.0'})
                            ylabel('Cortical depth (mm)')
                        elseif nCortical_Voxel == 45
                            Ylim = get(gca, 'ylim');
                            set(gca, 'YTick', linspace(Ylim(1), Ylim(end), 10));
                            yticklabels({'0','','0.5','','1.0','','1.5','','2.0','2.25'})
                            ylabel('Cortical depth (mm)')
                        else
                            Ylim = get(gca, 'ylim');
                            set(gca, 'YTick', linspace(Ylim(1), Ylim(end), 4));
                            yticklabels({'0','0.5','1.0','1.5','2.0','2.5','3.0'})
                            ylabel('Cortical depth (mm)')
                        end
                        
                    elseif Spatial_Resolution == 0.1
                        if nCortical_Voxel == 30
                            Ylim = get(gca, 'ylim');
                            set(gca, 'YTick', linspace(Ylim(1), Ylim(end), 6));
                        else
                            Ylim = get(gca, 'ylim');
                            set(gca, 'YTick', linspace(Ylim(1), Ylim(end), 4));
                        end
                        yticklabels({'0','0.5','1.0','1.5','2.0','2.5','3.0'})
                        ylabel('Cortical depth (mm)')
                    end
                    %                 set(gcf, 'Position', [500, 500, 1300, 1400]) %for ISMRM
                    %                 set(gcf, 'Position', [500, 500, 1300, 1400]) %for ISMRM
                    
                    set(gcf,'color','w');
                    set(gca,'FontSize',font_size,'FontWeight','Bold');
                    %                 set(gca,'FontSize',38,'FontWeight','Bold');
                    set(gca,'linewidth',axs_line_width);
                    saveas(h,strName105,'png');
                    
                    %%
                    %Detecting distinct fMRI profile, TR = 100ms
                    h=figure();
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
                    %     colormap=(sort(myColorOrder, 'descend'));
                    if (nCortical_Voxel == 60 && Spatial_Resolution == 0.05) || (nCortical_Voxel == 30 && Spatial_Resolution == 0.1)
                        c =colorbar('YTickLabel',{'0', '3'},'YTick', [0 1], 'Ydir','reverse');
                    elseif (nCortical_Voxel == 45 && Spatial_Resolution == 0.05)
                        c =colorbar('YTickLabel',{'0', '2.25'},'YTick', [0 1], 'Ydir','reverse');
                    else
                        c =colorbar('YTickLabel',{'0', '2'},'YTick', [0 1], 'Ydir','reverse');
                    end
                    c.Label.String = 'Cortical depth (mm)';
                    c.Label.Rotation = 90; % to rotate the text
                    
                    set(c,'YTickMode','manual');
                    ylabel('Percentage (%)');
                    if Filter_On == 0
                        max_value = 0.3;
                        %                     max_value = 0.15;
                        %                     min_value = -0.1;
                        %                     max_value = 0.1;
                        min_value = -0.05;
                        
                        ylim([min_value max_value])
                        yticklabels('manual');
                        yticks(min_value:0.05:max_value)
                        %                     yticklabels({'-10','-5','0','5','10','15','20','25','30'})
                        yticklabels({'-5','0','5','10','15','20','25','30'})
                    else
                        yticklabels('auto')
                    end
                    
                    xticklabels('manual')
                    if TR == 0.05
                        xlim([0 400])
                        xticks(0:100:400)
                    else
                        xlim([0 200])
                        xticks(0:50:200)
                    end
                    xticklabels({'0 sec','5 sec','10 sec','15 sec','20 sec'})
                    %                 xlabel('Time(sec)');
                    
                    strName1 = sprintf('Averaged Time Course');
                    
                    title(strName1);
                    saveas(h,strName1,'png');
                    
                    %%
                    % Get Standard Deviation of the baseline noise and
                    % Thresold signal with 2x SD Noise
                    noise_area_starting_point = 1;
                    noise_area_end_point = 10;
                    Threshold_2x_std_noise = 2*std(std(Threshold_SI_percentage(:,noise_area_starting_point:noise_area_end_point)));
                    Noise_Thresh_SI_percentage  = Threshold_SI_percentage .* (abs(Threshold_SI_percentage) > Threshold_2x_std_noise); %Hard Thresholding
                    
                    clear Onset_window_data;
                    Onset_window_data = Noise_Thresh_SI_percentage(:,11:30); % Thresholding, 0 to 2 mm
                    
                    [X3,Y3] = meshgrid(1:size(Onset_window_data,2), 1:size(Onset_window_data,1));
                    
                    %// Define a finer grid of points
                    [X4,Y4] = meshgrid(1:0.1:size(Onset_window_data,2), 1:0.1:size(Onset_window_data,1));
                    %// Interpolate the data and show the output
                    outData1 = interp2(X3, Y3, Onset_window_data, X4, Y4, 'Linear');
                    
                    h=figure();
                    imagesc((outData1),[0 0.02]);
                    
                    colormap(jet);
                    colorbar;
                    if TR == 0.05
                        xticks(1:100:400)
                    else
                        xticks(1:50:200);
                    end
                    xticklabels({'0s','0.5s','1.0s','1.5s','2.0s'});
                    %     xlabel('Time (sec)');
                    if Spatial_Resolution == 0.05
                        yticks(1:size(outData1,1)/6:size(outData1,1));
                    elseif Spatial_Resolution == 0.1
                        yticks(1:size(outData1,1)/6:size(outData1,1));
                    end
                    yticklabels({'0','0.5','1.0','1.5','2.0','2.5','3.0'});
                    ylabel('Cortical depth (mm)');
                    set(gcf,'color','w');
                    strName2 = sprintf('FP-S1 map');
                    title(strName2);
                    set(gca, 'XAxisLocation', 'top')
                    set(gca,'FontSize',font_size,'FontWeight','Bold');
                    set(gca,'linewidth',axs_line_width);
                    saveas(h,strName2,'png');
                    
                    %%
                    Temporal_SNR_mean = tSNR_function(abs(cortical_depth_map(:,:)));
                    if Filter_On == 0
                        save Temporal_SNR_mean Temporal_SNR_mean;
                    end
                    %%
                    %tSNR calculation
                    %tSNR calculation
                    %      for x = 1: size(total_cortical_depth_map,1)
                    %      start_idx = half_max_index - 10;
                    start_idx = half_max_index - 0.5/Spatial_Resolution;
                    if start_idx < 1
                        start_idx =1;
                    end
                    plot_linewidth = 2;
                    %                 for x = start_idx : start_idx + nCortical_Voxel +0.5/Spatial_Resolution -1
                    for x = start_idx : start_idx + 64 -1
                        if  x <= size(total_cortical_depth_map)
                            Temporal_SNR(x - start_idx + 1) = tSNR_function(abs(total_cortical_depth_map(x,:)));
                        end
                    end
                    
                    max_v = max(Temporal_SNR)
                    min_v = min(Temporal_SNR)
                    h=figure();
                    if num_image == 1
                        if nImage == num_image
                            plot(1:size(Temporal_SNR,2), Temporal_SNR, '-k','LineWidth',plot_linewidth);
                        else
                            plot(1:size(Temporal_SNR,2), Temporal_SNR, '-b','LineWidth',plot_linewidth);
                        end
                    elseif num_image == 2
                        plot(1:size(Temporal_SNR,2), Temporal_SNR, '-b','LineWidth',plot_linewidth);
                    else
                        plot(1:size(Temporal_SNR,2), Temporal_SNR, '-b','LineWidth',plot_linewidth);
                    end
                    set(gcf,'color','w');
                    strName1 = sprintf('tSNR plot');
                    title (strName1);
                    yticklabels('auto')
                    %      yticks(0:1:9);
                    if RS == 1 || TR == 1
                        ylim([0 max(Temporal_SNR)*1.1]);
                    else
                        ylim([0 50]);
                    end
                    %      yticklabels({'0','1','2','3', '4', '5', '6', '7', '8', '9'});
                    ylabel('Temporal SNR Value');
                    xlim([1 size(Temporal_SNR,2)])
                    %                 xlim([1 50])
                    %      xticks(0:5:size(total_cortical_depth_map,1))
                    %      array_index = 0:10:60;
                    array_index00 = 0:ceil(0.5/Spatial_Resolution):64;
                    array_index00(1) = 1;
                    xticks(array_index00)
                    if size(total_cortical_depth_map,1) == 64
                        xticklabels({'0','0.5','1.0','1.5','2.0','2.5','3.0'})
                    else
                        xticklabels({'0','0.5','1.0','1.5','2.0','2.5','3.0','3.5','4.0','4.5','5.0','5.5','6.0'})
                    end
                    %      xlim([ size(total_cortical_depth_map,1)]);
                    %      xticklabels({'0', '#1','#2','#3', '#4'});
                    %      xlabel('Voxel location from cortical surface (mm)');
                    xlabel('Cortical Depth (mm)');
                    %      legend('tSNR');
                    set(gca,'FontSize',font_size,'FontWeight','Bold');
                    set(gca,'linewidth',axs_line_width);
                    saveas(h,strName1,'png');
                    Temporal_SNR_plot_mean = mean(Temporal_SNR);
                    if Filter_On == 0
                        save Temporal_SNR Temporal_SNR;
                        save Temporal_SNR_plot_mean Temporal_SNR_plot_mean;
                    end
                    %%
                    
                    norm_cortical_map = zscore(squeeze(abs(cortical_depth_map(:,:))),0,2)';
                    
                    R = corrcoef(norm_cortical_map);
                    
                    Auto_Corr_Matrix_Setzero = R-diag(diag(R));
                    corr_min = -0.2;
                    corr_max = 1.0;
                    h=figure();
                    imagesc(Auto_Corr_Matrix_Setzero);
                    colormap(jet);
                    c = colorbar('Location','eastoutside','Ticks',[corr_min corr_max]);
                    c.Label.String = 'Correlation Coefficient';
                    c.Label.Rotation = 270; % to rotate the text
                    strName1 = sprintf('Voxel-by-voxel Corrleation');
                    
                    set(gca, 'XAxisLocation', 'top');
                    set(gcf,'color','w');
                    if nCortical_Voxel == 45
                        array_index3 = 0:ceil(0.25/Spatial_Resolution):nCortical_Voxel;
                        array_index3(1) = 1;
                        yticks(array_index3)
                        yticklabels({'0','','0.5','','1.0','','1.5','','2.0','2.25'})
                        ylabel('Cortical depth (mm)')
                        xticks(array_index3)
                        xticklabels({'0','','0.5','','1.0','','1.5','','2.0','2.25'})
                        xlabel('Cortical depth (mm)')
                    else
                        array_index3 = 0:ceil(0.5/Spatial_Resolution):nCortical_Voxel;
                        array_index3(1) = 1;
                        yticks(array_index3)
                        yticklabels({'0','0.5','1.0','1.5','2.0','2.5','3.0'})
                        ylabel('Cortical depth (mm)')
                        xticks(array_index3)
                        xticklabels({'0','0.5','1.0','1.5','2.0','2.5','3.0'})
                        xlabel('Cortical depth (mm)')
                    end
                    %                 caxis([-0.2 1.0]);% 10152018
                    caxis([corr_min corr_max]);
                    
                    axis image;
                    set(gca,'FontSize',font_size,'FontWeight','Bold');
                    set(gca,'linewidth',axs_line_width);
                    saveas(h,strName1,'png');
                    if Filter_On == 1
                        data_str1 = sprintf('Slice_%d_auto_correlation_matrix_',num_image);
                        saved_data_str1 = strcat(data_str1,data_num2str);
                        save (saved_data_str1,'Auto_Corr_Matrix_Setzero')
                    end
                    if Spatial_Resolution == 0.05
                        Average_voxel_layers(1,:)  = mean(cortical_depth_map(1:3,:),1);%L1,0~0.15 mm, averaged voxel
                        Average_voxel_layers(2,:)  = mean(cortical_depth_map(4:11,:),1);%L2/3,0.15~0.55 mm, averaged voxel
                        Average_voxel_layers(3,:)  = mean(cortical_depth_map(12:16,:),1);%L4,0.55~0.8 mm, averaged voxel
                        Average_voxel_layers(4,:)  = mean(cortical_depth_map(17:27,:),1);%L5,0.8~1.4 mm, averaged voxel
                        Average_voxel_layers(5,:)  = mean(cortical_depth_map(28:40,:),1);%L6,1.4~2.0 mm, averaged voxel
                        if nCortical_Voxel > 40
                            Average_voxel_layers(6,:)  = mean(cortical_depth_map(41:nCortical_Voxel,:),1);%WM,2.1~3.0 mm, averaged voxel
                        end
                    elseif Spatial_Resolution == 0.1
                        Average_voxel_layers(1,:)  = mean(cortical_depth_map(1:2,:),1);%L1,0~0.15 mm, averaged voxel
                        Average_voxel_layers(2,:)  = mean(cortical_depth_map(3:5,:),1);%L2/3,0.15~0.55 mm, averaged voxel
                        Average_voxel_layers(3,:)  = mean(cortical_depth_map(6:8,:),1);%L4,0.55~0.8 mm, averaged voxel
                        Average_voxel_layers(4,:)  = mean(cortical_depth_map(9:14,:),1);%L5,0.8~1.4 mm, averaged voxel
                        Average_voxel_layers(5,:)  = mean(cortical_depth_map(15:20,:),1);%L6,1.4~2.0 mm, averaged voxel
                        if nCortical_Voxel > 20
                            Average_voxel_layers(6,:)  = mean(cortical_depth_map(21:nCortical_Voxel,:),1);%WM,2.1~3.0 mm, averaged voxel
                        end
                    end
                    %% layer_by_lyaer signal
                    layer_plot_linewidth = 2;
                    %    plot(1:size(cortical_depth_map,2),zscore(mean(abs(cortical_depth_map(:,:)),1)),'-k','LineWidth',plot_linewidth);
                    h=figure();
                    for layer_idx = 1 : size(Average_voxel_layers,1)
                        subplot( size(Average_voxel_layers,1),1,layer_idx)
                        layer_signal(layer_idx,:) = zscore(abs(Average_voxel_layers(layer_idx,:)),0,2);
                        %                 layer_signal = (abs(Average_voxel_layers(layer_idx,:)));
                        plot(layer_signal(layer_idx,:),'-k','LineWidth',layer_plot_linewidth)
                        if layer_idx == 1
                            title('Layer #1');
                        elseif layer_idx == 2
                            title('Layer #2-3');
                        elseif layer_idx > 2 && layer_idx < 6
                            title(sprintf('Layer #%d',layer_idx+1));
                        elseif layer_idx == 6
                            title('WM');
                        end
                        set(gca,'FontSize',14,'FontWeight','Bold');
                        set(gca,'linewidth',axs_line_width);
                        xlim([0 size(Average_voxel_layers,2)]);
                        xticks(0:160*(1/TR):640*(1/TR))
                        xticklabels({'0 sec','160 sec','320 sec','480s','640 sec'})
                        ylim([-4 4])
                        
                    end
                    strName100='Layer-by-layer signal plot';
                    % title(strName);
                    set(gcf,'color','w');
                    set(gcf, 'Position', [500, 500, 800, 700])
                    saveas(h,strName100,'png');
                    %% layer_by_lyaer psd
                    layer_plot_linewidth = 2;
                    %    plot(1:size(cortical_depth_map,2),zscore(mean(abs(cortical_depth_map(:,:)),1)),'-k','LineWidth',plot_linewidth);
                    h=figure();
                    for layer_idx = 1 : size(Average_voxel_layers,1)
                        subplot(size(Average_voxel_layers,1),1,layer_idx)
                        %                 layer_signal(layer_idx,:) = zscore(abs(Average_voxel_layers(layer_idx,:)),0,2);
                        %                 layer_signal = (abs(Average_voxel_layers(layer_idx,:)));
                        [px_layer_signal,f] = pwelch(layer_signal(layer_idx,:),[],[],[],Fs, 'onesided');
                        %                 norm_layer_pxx = 10*log10(px_layer_signal)./ max(10*log10(px_layer_signal));
                        %                 norm_layer_pxx = px_layer_signal./ max(px_layer_signal);
                        plot(f,10*log10(px_layer_signal),'k','LineWidth',layer_plot_linewidth)
                        %                  plot(f,(norm_layer_pxx),'k','LineWidth',layer_plot_linewidth)
                        %                 plot(layer_signal,'-k','LineWidth',layer_plot_linewidth)
                        if layer_idx == 1
                            title('Layer #1');
                        elseif layer_idx == 2
                            title('Layer #2-3');
                        elseif layer_idx > 2 && layer_idx < 6
                            title(sprintf('Layer #%d',layer_idx+1));
                        elseif layer_idx == 6
                            title('WM');
                        end
                        set(gca,'FontSize',14,'FontWeight','Bold');
                        set(gca,'linewidth',axs_line_width);
                        %                 xlim([0 size(Average_voxel_layers,2)]);
                        %                 ylabel('Normalized PSD (a.u)');
                        %                 ylim([0 30])
                        %                 ylim([-0.2 1.2])
                        %     xticklabels('auto');
                        xticklabels('manual')
                        ylabel('PSD');
                        if Fs == 10
                            xlim([0 5])
                            xticks(0:1:5)
                            xticklabels({'0 Hz','1 Hz','2 Hz','3 Hz','4 Hz','5 Hz'})
                        elseif Fs == 20
                            xlim([0 10])
                            xticks(0:1:10)
                            xticklabels({'0 Hz','1 Hz','2 Hz','3 Hz','4 Hz','5 Hz','6 Hz','7 Hz','8 Hz','9 Hz','10 Hz'})
                        end
                        %                 if Fs == 10
                        %                     xlim([0.01 0.1])
                        %                     array_index8 = 0.01:0.01:0.1;
                        %                     %         array_index8(1) = 0.01;
                        %                     xticks(array_index8)
                        %                     xticklabels({'0.01','0.02','0.03','0.04','0.05', ...
                        %                         '0.06','0.07','0.08','0.09','0.1'})
                        %                 elseif Fs == 20
                        %                     xlim([0.01 0.1])
                        %                     xticks(0.01:0.01:0.1)
                        %                     xticklabels({'0.01','0.02','0.03','0.04','0.05', ...
                        %                         '0.06','0.07','0.08','0.09','0.1'})
                        %                 end
                        %                 xlabel('Frequency (Hz)');
                        %                 set(gcf, 'Position', [400, 400, 900, 350])
                        %                 set(gca,'FontSize',font_size,'FontWeight','Bold');
                        %                 set(gca,'linewidth',axs_line_width);
                        %                 set(gcf,'color','w');
                        %                 strName2 ='Zoomed & averaged PSD with whole time course';
                        %                 title(strName2);
                        %                 saveas(h,strName2,'png');
                    end
                    strName100='Layer-by-layer PSD plot';
                    % title(strName);
                    set(gcf,'color','w');
                    %             set(gcf, 'Position', [500, 500, 500, 800])
                    set(gcf, 'Position', [500, 500, 800, 700])
                    saveas(h,strName100,'png');
                    
                    %% layer_by_lyaer Zoomed PSD
                    layer_plot_linewidth = 2;
                    %    plot(1:size(cortical_depth_map,2),zscore(mean(abs(cortical_depth_map(:,:)),1)),'-k','LineWidth',plot_linewidth);
                    h=figure();
                    for layer_idx = 1 : size(Average_voxel_layers,1)
                        subplot(size(Average_voxel_layers,1) ,1,layer_idx)
                        %                 layer_signal = zscore(abs(Average_voxel_layers(layer_idx,:)),0,2);
                        %                 layer_signal = (abs(Average_voxel_layers(layer_idx,:)));
                        [px_layer_signal,f] = pwelch(layer_signal(layer_idx,:),[],[],[],Fs, 'onesided');
                        norm_layer_pxx = 10*log10(px_layer_signal)./ max(10*log10(px_layer_signal));
                        %                 norm_layer_pxx = px_layer_signal./ max(px_layer_signal);
                        %                  plot(f,10*log10(px_layer_signal),'k','LineWidth',layer_plot_linewidth)
                        plot(f,norm_layer_pxx,'k','LineWidth',layer_plot_linewidth)
                        %                 plot(layer_signal,'-k','LineWidth',layer_plot_linewidth)
                        if layer_idx == 1
                            title('Layer #1');
                        elseif layer_idx == 2
                            title('Layer #2-3');
                        elseif layer_idx > 2 && layer_idx < 6
                            title(sprintf('Layer #%d',layer_idx+1));
                        elseif layer_idx == 6
                            title('WM');
                        end
                        set(gca,'FontSize',14,'FontWeight','Bold');
                        set(gca,'linewidth',axs_line_width);
                        %                 xlim([0 size(Average_voxel_layers,2)]);
                        %                 ylabel('Normalized PSD (a.u)');
                        %                 ylim([0 30])
                        %                 ylim([-0.2 1.2])
                        %     xticklabels('auto');
                        xticklabels('manual')
                        if Fs == 10
                            xlim([0.01 0.5])
                            array_index9 = 0:0.1:0.5;
                            array_index9(1) = 0.01;
                            xticks(array_index9)
                            %          xticklabels({'0.01','0.1','0.2','0.3','0.4','0.5'})
                            xticklabels({'0.01','0.1','0.2','0.3','0.4','0.5'})
                        elseif Fs == 20
                            xlim([0.01 0.5])
                            array_index9 = 0:0.1:0.5;
                            array_index9(1) = 0.01;
                            xticklabels({'0.01','0.1','0.2','0.3','0.4','0.5'})
                        end
                        ylabel('Norm PSD');
                        %                 if Fs == 10
                        %                     xlim([0. 0.1])
                        %                     array_index8 = 0:0.01:0.1;
                        %                             array_index8(1) = 0.01;
                        %                     xticks(array_index8)
                        %                     xticklabels({'0.01','0.02','0.03','0.04','0.05', ...
                        %                         '0.06','0.07','0.08','0.09','0.1'})
                        %                 elseif Fs == 20
                        %                     xlim([0 0.1])
                        %                     xticks(0:0.01:0.1)
                        %                     xticklabels({'0.01','0.02','0.03','0.04','0.05', ...
                        %                         '0.06','0.07','0.08','0.09','0.1'})
                        %                 end
                        %                 xlabel('Frequency (Hz)');
                        %                 set(gcf, 'Position', [400, 400, 900, 350])
                        %                 set(gca,'FontSize',font_size,'FontWeight','Bold');
                        %                 set(gca,'linewidth',axs_line_width);
                        %                 set(gcf,'color','w');
                        %                 strName2 ='Zoomed & averaged PSD with whole time course';
                        %                 title(strName2);
                        %                 saveas(h,strName2,'png');
                    end
                    strName101='Zoomed Layer-by-layer PSD (0p01-0p1) plot';
                    % title(strName);
                    set(gcf,'color','w');
                    %             set(gcf, 'Position', [500, 500, 500, 800])
                    set(gcf, 'Position', [500, 500, 600, 700])
                    saveas(h,strName101,'png');
                end%%this for loop is for image number (e.g. 1 or 3)
            end %%for filtering case
        end % for exp_case_num
    end %resting or task
end% for animal number
%for next job
% match_acq_fmri_corr_ms_2slice_bpf_compen_coh_Auto_20181026_vF
