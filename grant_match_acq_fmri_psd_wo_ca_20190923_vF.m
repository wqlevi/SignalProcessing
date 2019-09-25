clear;
close all;
% clc;


addpath('/Users/Daniel/2_Bruker_Scanner/Calcium_by_Filip');
addpath('/Users/Daniel/2_Bruker_Scanner/utils');
% path = '/Users/Daniel/2_Bruker_Scanner/1_Bruker_Data/Biopac/02232019_line/Matfile';
% path = '/Users/Daniel/2_Bruker_Scanner/1_Bruker_Data/2018yrs/03122018/Biopac_Matfile/';
% path = '/Users/Daniel/2_Bruker_Scanner/1_Bruker_Data/2018yrs/05062018/Biopac_Matfile/';
% path ='/Users/Daniel/2_Bruker_Scanner/1_Bruker_Data/04062019_linescan/Biopac/linescan_TK/';
% path = '/Users/Daniel/2_Bruker_Scanner/1_Bruker_Data/2018yrs/09062018/Biopac_Matfile/';
path = '/Users/Daniel/2_Bruker_Scanner/1_Bruker_Data/2018yrs/09252018/Biopac_Matfile/';
% path = '/Users/Daniel/2_Bruker_Scanner/1_Bruker_Data/2018yrs/08162018/Biopac_Matfile/';
% linedata_path ='/Users/Daniel/2_Bruker_Scanner/1_Bruker_Data/2018yrs/09062018/1.Pe1/';
% linedata_path ='/Users/Daniel/2_Bruker_Scanner/1_Bruker_Data/2018yrs/08162018/1.OT1/';
% linedata_path = '/Users/Daniel/2_Bruker_Scanner/1_Bruker_Data/2018yrs/09252018/1.Px1/';
path_root = '/Users/Daniel/2_Bruker_Scanner/1_Bruker_Data/2018yrs/SS_LINE/';
linedata_path = strcat(path_root,'/09252018/');
% path = '/Users/Daniel/2_Bruker_Scanner/1_Bruker_Data/2018yrs/10062018/Biopac_Matfile/';
% linedata_path ='/Users/Daniel/2_Bruker_Scanner/1_Bruker_Data/2018yrs/10062018/1.PI1/';

font_size = 32;
axs_font_size = 40;%ppt font 14
set(0,'defaultAxesFontName', 'Arial')
set(0,'defaultTextFontName', 'Arial')
set(0,'defaultTextFontSize', font_size)
set(0,'defaultAxesFontSize',axs_font_size)
set(0,'defaultTextFontWeight', 'Bold')
set(0,'defaultAxesFontWeight', 'Bold')

axs_line_width = 2;
plot_linewidth = 2;
n_plots = 50;

SAVE_On = 0;
RegressPoly_On = 0; %it doesn't work with evoked fMRI
RegressPhysio_On = 0;
card_f_low = 6;
card_f_high = 7;
Cardiac_Reg_On = 0;
Respiration_Reg_On = 0;

% Ca_On = 0; %1 :Calcium, 0 : without Calcium

Excel_Export_On = 0;
Filter_On = 1; %calculate specific frequency power, 0.01-0.1 Hz

TR50 = 1;

if TR50 == 1
    TR = 0.05; %s
    b4trig = 20; %prestimbaseline
    nSlice = 400;
%     total_scan_time = 640; %second
else
    TR = 0.1; %s
    b4trig = 10; %prestimbaseline
    nSlice = 200;
%     total_scan_time = 640; %second
end
% Total_Slice = 3;
Total_Slice = 1;
Spatial_Resolution = 0.05;
% TR = 0.1; %
% TR = 0.05; %
% nSlice = 200;
% TR = 0.06; %

nY = 32;
TA = nSlice*TR*nY;

% b4trig = 10; %prestimbaseline
% b4trig = 20; %prestimbaseline
total_scan_time = TA; %second
fmri_dummy = ones(total_scan_time*(1/TR), 1);


Fs = floor(1/TR);
if Fs <  20
    Fs1= floor(1/(TR/2)); %for check higher freq than nyquist
% fmri_dummy = ones(total_scan_time*(1/TR),1);
else 
    Fs1= Fs;
end

% Spatial_Resolution = 0.1; %0.1
Spatial_Resolution = 0.05; %0.05
if Spatial_Resolution == 0.1
     nCortical_Voxel = 20;
elseif Spatial_Resolution == 0.05
    nCortical_Voxel = 40; 
%     nCortical_Voxel = 45; 
end

RS = 0;
% Exp_Number_Array1 = [24 29 31 35 37 39 41 44 46 48 50 52 54 56 58 60 62 64 66 68 70 72 74]; %Task-related fMRI
% Exp_Number_Array1 = [72 74]; %Task-related fMRI
% Exp_Number_Array1 = [41, 43, 48, 51]; %Rest-state fMRI
% Exp_Number_Array1 = [38]; %evoked fMRI
% Exp_Number_Array1 = [29]; %evokedfMRI, TR 100 ms

Exp_Number_Array1 = [31]; %evokedfMRI, TR 50 ms,09252018

% Exp_Number_Array1 = [33]; %rest-state fMRI, TR 50 ms


for exp_case_num = 1 : size(Exp_Number_Array1,2)

    data_num2str = num2str(Exp_Number_Array1(exp_case_num));
    data_num = Exp_Number_Array1(exp_case_num);

    %load cortical depth map data
    temp_data_path = strcat(linedata_path,'/',data_num2str,'/');
%     cd(temp_data_path);
    %% filtering fMRI signals
    clear linescanning_data;
    if Filter_On == 0 %without filtering, original signal
        for num_slice = 1 : Total_Slice
            if num_slice == 1
                strTemp = 'Results_Slice1/';
                cd(strcat(temp_data_path,strTemp));
                linescanning_data(:,:,num_slice) = cell2mat(struct2cell(load('cortical_depth_map_s1.mat'))); %load data
            elseif num_slice == 2
                strTemp = 'Results_Slice2/';
                cd(strcat(temp_data_path,strTemp));
                linescanning_data(:,:,num_slice) = cell2mat(struct2cell(load('cortical_depth_map_s2.mat'))); %load data
            elseif num_slice == 3
                strTemp = 'Results_Slice3/';
                cd(strcat(temp_data_path,strTemp));
                linescanning_data(:,:,num_slice) = cell2mat(struct2cell(load('cortical_depth_map_s3.mat'))); %load data
            end
            cortex{num_slice}{exp_case_num} = abs(squeeze(linescanning_data(:,:,num_slice)))'; 
        end
    elseif Filter_On == 1 %band pass filtering, 0.01-0.1 Hz
        strTemp2 = 'Results_Filtering_0p1';
        for num_slice = 1 : Total_Slice
            if num_slice == 1
                strTemp = 'Results_Slice1/';
                cd(strcat(temp_data_path,strTemp,strTemp2));
                linescanning_data(:,:,num_slice) = cell2mat(struct2cell(load('cortical_depth_map_s1.mat'))); %load data
            elseif num_slice == 2
                strTemp = 'Results_Slice2/';
                cd(strcat(temp_data_path,strTemp,strTemp2));
                linescanning_data(:,:,num_slice) = cell2mat(struct2cell(load('cortical_depth_map_s2.mat'))); %load data
            elseif num_slice == 3
                strTemp = 'Results_Slice3/';
                cd(strcat(temp_data_path,strTemp,strTemp2));
                linescanning_data(:,:,num_slice) = cell2mat(struct2cell(load('cortical_depth_map_s3.mat'))); %load data
            end
            cortex{num_slice}{exp_case_num} = abs(squeeze(linescanning_data(1:nCortical_Voxel,:,num_slice)))'; 
        end
    end
    %% replace cortical voxels
    clear linescanning_data;
    for num_slice = 1 : Total_Slice
        linescanning_data(1:nCortical_Voxel,:,num_slice) = cortex{num_slice}{exp_case_num}(:,1:nCortical_Voxel)' ;
    end
    
      
    %%
    clear biopac_data;
    clear data2 beg fin;
    cd(path);
     
    if RS == 0
%         cd(strcat(path,'\0_Task_related'));
        %load biopac data
        biopac_data= load(strcat('scan_',data_num2str,'_TK.mat')); 
    else
%         cd(strcat(path,'\1_Rest_State'));
        biopac_data= load(strcat('scan_',data_num2str,'_RS.mat')); %load data
    end
    %         mkdir Results_Figures;
%     cd(path);
    cd(temp_data_path);
    folder_str = strcat('Results_Figures_Physio_#',data_num2str);
    mkdir(folder_str);
    pathstr = strcat(temp_data_path,'/', folder_str);
    cd(pathstr); %to save figures for cardiac
    if Filter_On == 1 && RegressPhysio_On == 0 %band pass filtering only, 0.01-0.1 Hz
        folder_str2 = strcat('Results_Filtering_0p1_#',data_num2str);
        mkdir(folder_str2);
        pathstr2 = strcat(pathstr,'/', folder_str2);
        cd(pathstr2); %to save figures.
    elseif Filter_On == 0 && RegressPhysio_On == 1
        folder_str2 = strcat('Results_RegressPhysio_#',data_num2str);
        mkdir(folder_str2);
        pathstr2 = strcat(pathstr,'/', folder_str2);
        cd(pathstr2); %to save figures.
    elseif Filter_On == 1 && RegressPhysio_On == 1
        folder_str2 = strcat('Results_Filter0p1_RegressPhysio_#',data_num2str);
        mkdir(folder_str2);
        pathstr2 = strcat(pathstr,'/', folder_str2);
        cd(pathstr2); %to save figures.
    end
%%      clear linescanning_data;
     for num_slice = 1 : Total_Slice 
         % just for checking the line-scanning signals
         Mean_linescanning_fMRI = zscore(mean(abs(linescanning_data(:,:,num_slice)),1),0,2);
         h=figure();
         if num_slice == 1
             if num_slice == Total_Slice
                 plot(1:size(Mean_linescanning_fMRI,2),Mean_linescanning_fMRI,'-k','LineWidth',plot_linewidth);
             else
                 plot(1:size(Mean_linescanning_fMRI,2),Mean_linescanning_fMRI,'-b','LineWidth',plot_linewidth);
             end
         elseif num_slice == 2
             plot(1:size(Mean_linescanning_fMRI,2),Mean_linescanning_fMRI,'-r','LineWidth',plot_linewidth);
         else
             plot(1:size(Mean_linescanning_fMRI,2),Mean_linescanning_fMRI,'-k','LineWidth',plot_linewidth);
         end
         strName1 = sprintf('Line-scanning fMRI, averaged voxel time course');
         %     title(strName1);
         xlim([0 640*(1/TR)])
         xticks(0:160*(1/TR):640*(1/TR))
         xticklabels({'0s','160s','320s','480s','640s'})
         set(gca,'xtick',[]);
         %     xlabel('Time(sec)');
         ylim([-4.5 4.5])
         ylabel('fMRI Signal (a.u)', 'FontSize', font_size );
         ax1 = gca;
         % ax.Visible = 'off'
         ax1.YLabel.Visible = 'on'; % remove y-axis
         %     ax1.XAxis.Visible = 'on'; % remove x-axis
         %    set(gca,'ylabel','');
         %     title(strName1,'FontSize', font_size );
         set(gcf,'color','w');
         %     set(gca,'FontWeight','Bold');
         set(gca,'linewidth',axs_line_width);
         set(gcf, 'Position', [500, 500, 900, 300])
         strName12 = sprintf('Slice #%d, Origin Averaged voxel time course',num_slice);
         saveas(h,strName12,'png'); %save line-scanning figure
         %%
         if Filter_On == 1
             cortical_depth_map = linescanning_data;
             norm_cortical_map_temp = abs(cortical_depth_map(:,:))./max(max(abs(cortical_depth_map(:,:))));
             Threshold_SI_MAP = abs(cortical_depth_map(:,:))./mean(norm_cortical_map_temp(:,:),2);
             norm_cortical_map = abs(Threshold_SI_MAP(:,:))./max(max(abs(Threshold_SI_MAP(:,:))));
             %                     norm_cortical_map = (abs(Threshold_SI_MAP(:,:))-min(min(abs(Threshold_SI_MAP(:,:)))))./(max(max(abs(Threshold_SI_MAP(:,:))))-min(min(abs(Threshold_SI_MAP(:,:)))));
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
                     caxis([0.3 1])
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
%              title(strName0);
             set(gcf,'color','w');
             set(gca,'FontSize',font_size,'FontWeight','Bold');
             set(gca,'linewidth',axs_line_width);
             
             xlim([0 TA*(1/TR)])
             if TA == 640
                 xticks(0:TA/4*(1/TR):TA*(1/TR))
                 xticklabels({'0s','160s','320s','480s','640s'})
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
                     set(gca, 'YTick', linspace(Ylim(1), Ylim(end), 5));
                     yticklabels({'0','0.5','1.0','1.5','2.0'})
%                      ylabel('Cortical depth (mm)')
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
             set(gcf, 'Position', [500, 500, 1000, 300])
             %             set(gcf, 'Position', [500, 500, 2500, 500])
             saveas(h,strName0,'png');
         end
         %% zoomed line
         zoomed_plot_linewidth =4;
         Zoomed_Mean_linescanning_fMRI = Mean_linescanning_fMRI(:,1:10*(1/TR));
         h=figure();
         if num_slice == 1
             if num_slice == Total_Slice
                 plot(1:size(Zoomed_Mean_linescanning_fMRI,2),Zoomed_Mean_linescanning_fMRI,'-k','LineWidth',zoomed_plot_linewidth);
             else
                 plot(1:size(Zoomed_Mean_linescanning_fMRI,2),Zoomed_Mean_linescanning_fMRI,'-b','LineWidth',zoomed_plot_linewidth);
             end
         elseif num_slice == 2
             plot(1:size(Zoomed_Mean_linescanning_fMRI,2),Zoomed_Mean_linescanning_fMRI,'-r','LineWidth',zoomed_plot_linewidth);
         else
             plot(1:size(Zoomed_Mean_linescanning_fMRI,2),Zoomed_Mean_linescanning_fMRI,'-k','LineWidth',zoomed_plot_linewidth);
         end
         strName1 = sprintf('Zoomed line-scanning fMRI');
         %     title(strName1);
         xlim([0 10*(1/TR)])
         xticks(0:1*(1/TR):10*(1/TR))
         xticklabels({'0s','1s','2s','3s','4s','5s','6s','7s','8s','9s','10s'})
         %     xlabel('Time(sec)');
         ylim([-4.5 4.5])
         ylabel('fMRI Signal (a.u)')
         ax1=gca;
         ax1.YLabel.Visible = 'off'; % remove y-axis
         %     title(strName1);
         set(gca,'xtick',[]);
         set(gcf,'color','w');
         %     set(gca,'FontSize',font_size,'FontWeight','Bold');
         set(gca,'linewidth',axs_line_width);
         set(gcf, 'Position', [500, 500, 900, 300])
         strName12 = sprintf('Slice #%d, Zoomed origin Averaged voxel time course',num_slice);
         saveas(h,strName12,'png'); %save line-scanning figur
     end
     
     %% extract Linescanning_fMRI Excel file
%      [num,txt,raw] = xlsread('Table2');
    if Excel_Export_On == 1
%         delete Linescanning_fMRI_TR50ms.xlsx
        filename = 'Linescanning_fMRI_TR50ms.xlsx';
        A1 = Mean_linescanning_fMRI';
        T1 =array2table(A1);
        T1.Properties.VariableNames = {'Mean_linescanning_fMRI'};

        B1 = Zoomed_Mean_linescanning_fMRI';
        T2 =array2table(B1);
        T2.Properties.VariableNames = {'Zoomed_Mean_linescanning_fMRI'};

        writetable(T1,filename,'Sheet', 1,'Range','A1') ;
        writetable(T2,filename,'Sheet', 1,'Range','B1') ;
    end
%     fid = fopen('Linescanning_fMRI_TR50ms.xlsx');
%     fgetl(fid);
%     fclose(fid);
%     writetable(T2,filename,'Sheet', 2) ;
    %% match_acq_fmri fucntion made by filip
    [data2,fmri_out,beg,fin] = match_acq_fmri(biopac_data,fmri_dummy,TR,b4trig);
%%
% your code here
    chan_card = 1; %for blood pressure
    chan_resp = 6; %for respiration
    chan_stim = 5; %for stimulation
    chan_ca = 4; %for stimulation
    cardiac_match{exp_case_num}       = (data2.channels{chan_card}.data).';
    respiratory_match{exp_case_num}    = (data2.channels{chan_resp}.data).';
    stim_match{exp_case_num}          = (data2.channels{chan_stim}.data).';%%%
    ca_match{exp_case_num}          = (data2.channels{chan_ca}.data).';%%%
    samples_per_second    = (data2.channels{chan_card}.samples_per_second).';%%% blood
    samples_per_second2   = (data2.channels{chan_resp}.samples_per_second).';%%%respir
    samples_per_second3   = (data2.channels{chan_stim}.samples_per_second).';%%%stim

    
    blood_press_temp      =  cardiac_match{exp_case_num};
    repiratory_press_temp =  respiratory_match{exp_case_num};
    stim_temp             =  stim_match{exp_case_num};%%%
    
%     samples_per_second    = (data2.channels{chan_card}.samples_per_second).';%%% blood
%     samples_per_second2   = (data2.channels{chan_resp}.samples_per_second).';%%%respir
    
     data_num = num2str(Exp_Number_Array1(exp_case_num));
    
    h=figure();
    subplot (3,1,1);
    plot(1:size(blood_press_temp,1),blood_press_temp(:,1),'r-');
    title('blood press')
    subplot (3,1,2);
    plot(1:size(repiratory_press_temp,1),repiratory_press_temp(:,1),'b-');
    title('repiratory freq')
    subplot (3,1,3);
    plot(1:size(stim_temp,1),stim_temp(:,1),'m-');
    title('Stimulation')
%     subplot (4,1,4);
%     plot(1:size(Ca,1),Ca(:,1),'k-');
%     title('Calcuim')
    if RS == 1
        strName3 = sprintf('Rest #%d, original biopac',str2num(data_num));
    else
        strName3 =  sprintf('Evoked #%d, original biopac',str2num(data_num));
    end
    set(gcf,'color','w');
    saveas(h,strName3,'png');
    
    %% Original Blood Pressure
    h=figure();
    res_blood_temp =  zscore(blood_press_temp,0,1);
    plot(1:size(res_blood_temp,1),res_blood_temp(:,1),'r-','LineWidth',plot_linewidth);
    % title('Blood Pressure')
    ylabel('Signal (a.u)');
    ax1=gca;
    ax1.YLabel.Visible = 'off'; % remove y-axis
    xlim([0 size(res_blood_temp,1)])
    xticklabels('manual')
    Xlim = get(gca, 'xlim');
    set(gca, 'XTick', linspace(Xlim(1), Xlim(end), 5));
    if TA == 640
        xticklabels({'0s','160s','320s','480s','640s'})
    elseif TA == 768
        xticklabels({'0s','192s','384s','576s','768s'})
    end
    xticklabels('manual')
    yticklabels('auto')
    ylim([-4 4])
    if RS == 1
        strName3 = sprintf('Rest #%d, original blood press biopac',str2num(data_num));
    else
        strName3 =  sprintf('Evoked #%d, original blood press biopac',str2num(data_num));
    end
    set(gcf,'color','w');
     set(gca,'xtick',[]);
    % set(gca,'FontSize',font_size,'FontWeight','Bold');
    set(gca,'linewidth',axs_line_width);
    set(gcf, 'Position', [500, 500, 1200, 300])
    saveas(h,strName3,'png');
    %% export blood pressure to excel
    if Excel_Export_On == 1
        filename = 'Blood_Pressure_TR50ms.xlsx';
        blood_press = single(reshape(res_blood_temp, [TA samples_per_second]));% total scan time 640 sec x sample per second 5000
        T =array2table(blood_press);

        %     T.Properties.VariableNames = {'Blood_Pressure'};
        writetable(T,filename,'Sheet', 1) ;
    end
    %% Zoomed Blood
    zoomed_plot_linewidth =4;
    Second = 10;%10sec
    nWindow = samples_per_second *Second;
    min_x = 0;
    max_x = nWindow;
    h=figure();
    zoomed_res_blood_temp =  zscore(blood_press_temp(1:nWindow),0,1);
    plot(1:size(zoomed_res_blood_temp,1),zoomed_res_blood_temp(:,1),'r-','LineWidth',zoomed_plot_linewidth);
    % title('Zoomed Blood Pressure')
   
    
    xlim([min_x  max_x ])
    % patch([min_x  max_x], [-4.5  -4.5 4.5 4.5], [0.5 0.8 0.9])
    
    xticklabels('manual')
    yticklabels('auto')
    ylabel('Signal (a.u)');
    ylim([-4.5 4.5])
    if RS == 1
        strName3 = sprintf('Rest #%d, zoomed original blood press biopac',str2num(data_num));
    else
        strName3 =  sprintf('Evoked #%d, zoomed original blood press biopac',str2num(data_num));
    end
    set(gcf,'color','w');
    % set(gca, 'color', [0.9 0.9 0.9])
    % set(gca,'FontSize',font_size,'FontWeight','Bold');
    set(gca,'linewidth',axs_line_width);
    % set(gca,'visible','off')
    % set(gca,'xtick',[])
    % set(gca,'XColor', 'none')
    ax1=gca;
    ax1.YLabel.Visible = 'off'; % remove y-axis
    % ax.Visible = 'off'
    ax1.YAxis.Visible = 'on'; % remove y-axis
    ax1.XAxis.Visible = 'on'; % remove x-axis
    set(gcf, 'Position', [500, 500, 900, 300])
    saveas(h,strName3,'png');
    
    %% export blood pressure to excel
    if Excel_Export_On == 1
        zoomed_blood_press = single(reshape(zoomed_res_blood_temp, [Second samples_per_second]));
        T1 =array2table(zoomed_blood_press);
        %     T.Properties.VariableNames = {'Blood_Pressure'};
        writetable(T1,filename,'Sheet', 2) ;
    end
%     xlswrite('blood_pressure.xlsx',res_blood_temp');
%     a=csvread('blood_pressure.csv');
%     temp_a = a';
%     b=regexp(temp_a,',','split');
%     c=reshape([b{:}],numel(b{1}),[])';
%     xlswrite('test.xlsx',c);
%% Original Respiration

    h=figure();
    res_repiratory_temp =  zscore(repiratory_press_temp,0,1);
    plot(1:size(res_repiratory_temp,1),res_repiratory_temp(:,1),'b-','LineWidth',plot_linewidth);
    % title('Respiration')
    xlim([0 size(res_repiratory_temp,1)])
    xticklabels('manual')
    Xlim = get(gca, 'xlim');
    set(gca, 'XTick', linspace(Xlim(1), Xlim(end), 5));
    if TA == 640
        xticklabels({'0s','160s','320s','480s','640s'})
    elseif TA == 768
        xticklabels({'0s','192s','384s','576s','768s'})
    end

    yticklabels('auto')
    ylabel('Signal (a.u)');
    ax1=gca;
    ax1.YLabel.Visible = 'off'; % remove y-axis
    ylim([-3 3])
    if RS == 1
        strName3 = sprintf('Rest #%d, original respiration biopac',str2num(data_num));
    else
        strName3 =  sprintf('Evoked #%d, original respiration biopac',str2num(data_num));
    end
    set(gcf,'color','w');
    % set(gca,'FontSize',font_size,'FontWeight','Bold');
    set(gca,'linewidth',axs_line_width);
    set(gcf, 'Position', [500, 500, 900, 300])
    saveas(h,strName3,'png');
    
     %% export respiration to excel
     if Excel_Export_On == 1
         filename = 'Respiration_TR50ms.xlsx';
         respiration = single(reshape(res_repiratory_temp, [TA samples_per_second2]));% total scan time 640 sec x sample per second 5000
         T =array2table(respiration);
         
         %     T.Properties.VariableNames = {'Blood_Pressure'};
         writetable(T,filename,'Sheet', 1) ;
     end
    %% Zoomed Respiration
    % zoomed_plot_linewidth =1.5;
    h=figure();
    
    nWindow = samples_per_second2 *Second;%10sec
    min_x = 0;
    max_x = nWindow;
    
    zoomed_res_repiratory_temp =  zscore(repiratory_press_temp(1:nWindow,1),0,1);
    plot(1:size(zoomed_res_repiratory_temp,1),zoomed_res_repiratory_temp(:,1),'b-','LineWidth',zoomed_plot_linewidth);
    % title('Zoomed Respiration')
    
    xlim([min_x  max_x ])
    % patch([min_x  max_x], [-4.5  -4.5 4.5 4.5], [0.5 0.8 0.9])

    xticklabels('manual')
    Xlim = get(gca, 'xlim');
    set(gca, 'XTick', linspace(Xlim(1), Xlim(end), 11));
    if TA == 640
    %     xticklabels({'0s','160s','320s','480s','640s'})
        xticklabels({'0s','1s','2s','3s','4s','5s','6s','7s','8s','9s','10s'})
    end

    yticklabels('auto')
    ylabel('Signal (a.u)');
    ylim([-3 3])
    if RS == 1
        strName3 = sprintf('Rest #%d, zoomed original respiration biopac',str2num(data_num));
    else
        strName3 =  sprintf('Evoked #%d, zoomed original respiration biopac',str2num(data_num));
    end
    set(gcf,'color','w');
    % set(gca,'FontSize',font_size,'FontWeight','Bold');
    set(gca,'linewidth',axs_line_width);
    % set(gca,'visible','off')
    % set(gca,'xtick',[])
    % set(gca,'XColor', 'none')
    % ax1 = gca;
    ax1=gca;
    ax1.YLabel.Visible = 'off'; % remove y-axis
    % ax.Visible = 'off'
    ax1.YAxis.Visible = 'on'; % remove y-axis
    ax1.XAxis.Visible = 'on'; % remove x-axis
    set(gcf, 'Position', [500, 500, 900, 300])
    saveas(h,strName3,'png');
    %%
    if Excel_Export_On == 1
        zoomed_respiration = single(reshape(zoomed_res_repiratory_temp, [Second samples_per_second2]));
        T1 =array2table(zoomed_respiration);
        %     T.Properties.VariableNames = {'Blood_Pressure'};
        writetable(T1,filename,'Sheet', 2) ;
    end
     %%
     %2)Downsampling
    synchorized_blood_press_temp =resample(blood_press_temp,Fs,samples_per_second)';
    synchorized_respiratory_press_temp =resample(repiratory_press_temp,Fs,samples_per_second2)';
    synchorized_stim_temp =resample(stim_temp,Fs,samples_per_second3)';
    
    %add offset tricky, 25 08 2018
%     synchorized_blood_press(1,1) = synchorized_blood_press(1,2);
%     synchorized_respiratory_press(1,1) = synchorized_respiratory_press(1,2);
    
   
    h=figure();
    nSingal = 2;
    for index = 1 : nSingal
        subplot (nSingal,1,index);
        if index == 1
            %for blood pressure
            res_blood =  zscore(synchorized_blood_press_temp,0,2);      
%             plot(1:size(synchorized_blood_press,2),synchorized_blood_press(1,:),'r-');
            plot(1:size(res_blood,2),res_blood(1,:),'r-');
            strName1 = sprintf('Synchronized Blood Pressure');
            title(strName1);
        elseif index == 2
            %for respiratory press
            res_respiratory =  zscore(synchorized_respiratory_press_temp,0,2);   
%             plot(1:size(synchorized_respiratory_press,2),synchorized_respiratory_press(1,:),'b-');
            plot(1:size(res_respiratory,2),res_respiratory(1,:),'b-');
            strName1 = sprintf('Synchronized Respiratory Signal');
            title(strName1);
%         elseif index == 3
%             %for stimulation signal
%             res_stim =  zscore(synchorized_stim_temp,0,2);   
% %             plot(1:size(synchorized_stim,2),synchorized_stim(1,:),'m-');
%             plot(1:size(res_stim,2),res_stim(1,:),'m-');
%             strName1 = sprintf('Synchronized Stimulation Signal');
%             title(strName1);
%             ylim([-0.07 0.07])
        end
        ylim([-5 5])
        yticklabels('auto')
        xticklabels('manual')
        xlim([0 TA*(1/TR)])
        if TA == 640
            xticks(0:TA/4*(1/TR):TA*(1/TR))
            xticklabels({'0s','160s','320s','480s','640s'})
        elseif TA == 768
            xticks(0:TA/4*(1/TR):TA*(1/TR))
            xticklabels({'0s','192s','384s','576s','768s'})
        end
        %     xlabel('Time(sec)');
        set(gcf,'color','w');
        set(gca,'FontSize',font_size,'FontWeight','Bold');
        set(gca,'linewidth',axs_line_width);
        set(gcf, 'Position', [500, 500, 900, 300])
        if RS == 1
            strName4 = sprintf ('Rest #%d, Synchronized signals', str2num(data_num));
        else
            strName4 = sprintf ('Evoked #%d, Synchronized signals', str2num(data_num));
        end
        saveas(h,strName4,'png');
    end
    
    %%
    res_line(:,:) = squeeze(mean(abs(linescanning_data(:,:)),1));
%             res_line2 = zscore(abs(squeeze(res_line(:,:,num_slice))),0,2);
    res_line2 = res_line(:,:) - mean(res_line(:,:),2);
%     x =res_line2;
    [pxx,f] = pwelch(res_line2,[],[],[],Fs, 'onesided');
    h= figure();
    plot(f,10*log10(pxx),'k','LineWidth',plot_linewidth)
    ylabel('PSD');
    ax1=gca;
    ax1.YLabel.Visible = 'off'; % remove y-axis
    ylim([-30 20]);
    
    xticklabels('manual')
    Xlim = get(gca, 'xlim');
    if Fs == 10
        set(gca, 'XTick', linspace(Xlim(1), Xlim(end), 6));
        xticklabels({'0','1','2','3','4','5'})
    elseif Fs == 20
        set(gca, 'XTick', linspace(Xlim(1), Xlim(end),11));
        xticklabels({'0','1','2','3','4','5','6','7','8','9','10'})
    end
    ax1.XLabel.Visible = 'off'; % remove y-axis
    set(gca,'xtick',[]);
%     yticklabels('auto')
%     ylabel('Signal (a.u)');
%     ax1=gca;
%     ax1.YLabel.Visible = 'off'; % remove y-axis
    xlabel('Frequency (Hz)');
    set(gcf, 'Position', [500, 500, 900, 300])
    set(gca,'FontSize',font_size,'FontWeight','Bold');
    set(gca,'linewidth',axs_line_width);
    set(gcf,'color','w');
    strName2 = sprintf('Slice #%d, Averaged Line-scanning PSD',num_slice);
%     title(strName2);
    saveas(h,strName2,'png');
    %% Blood PSD
    h=figure();
    %for blood pressure
    x1 = synchorized_blood_press_temp(1,:);
    [pxx_blood,f] = pwelch(x1,[],[],[],Fs1, 'onesided');
    plot(f,10*log10(pxx_blood),'r-','LineWidth',plot_linewidth);
    strName1 = sprintf('Blood Pressure PSD');
    %      title(strName1);
    
    ylabel('PSD');
    ax1=gca;
    ax1.YLabel.Visible = 'off'; % remove y-axis
    ylim([-60 60]);
    
    xticklabels('manual')
    Xlim = get(gca, 'xlim');
    if Fs == 10
        set(gca, 'XTick', linspace(Xlim(1), Xlim(end), 6));
        xticklabels({'0','1','2','3','4','5'})
    elseif Fs == 20
        set(gca, 'XTick', linspace(Xlim(1), Xlim(end),11));
        xticklabels({'0','1','2','3','4','5','6','7','8','9','10'})
    end
    ax1.XLabel.Visible = 'off'; % remove y-axis
    set(gca,'xtick',[]);

     set(gcf,'color','w');
     set(gca,'FontSize',font_size,'FontWeight','Bold');
     set(gca,'linewidth',axs_line_width);
     set(gcf, 'Position', [500, 500, 700, 300])
     
     strName4 = sprintf('Blood Pressure #%d, PSD plots',str2num(data_num));
     saveas(h,strName4,'png');

    %%
    h=figure();
    
    %for respiratory press
    x2 = synchorized_respiratory_press_temp(1,:);
    [pxx_resp,f] = pwelch(x2,[],[],[],Fs1, 'onesided');
    plot(f,10*log10(pxx_resp),'b-','LineWidth',plot_linewidth);
    strName1 = sprintf('Respiratory Signal PSD');
%     title(strName1);
    ylabel('PSD');
    ax1=gca;
    ax1.YLabel.Visible = 'off'; % remove y-axis
    ylim([-60 60]);
    
    xticklabels('manual')
    Xlim = get(gca, 'xlim');
    if Fs == 10
        set(gca, 'XTick', linspace(Xlim(1), Xlim(end), 6));
        xticklabels({'0','1','2','3','4','5'})
    elseif Fs == 20
        set(gca, 'XTick', linspace(Xlim(1), Xlim(end),11));
        xticklabels({'0','1','2','3','4','5','6','7','8','9','10'})
    end
    ax1.XLabel.Visible = 'off'; % remove y-axis
%     set(gca,'xtick',[]);
    xlabel('Frequency (Hz)');
    set(gcf,'color','w');
    set(gca,'FontSize',font_size,'FontWeight','Bold');
    set(gca,'linewidth',axs_line_width);
    set(gcf, 'Position', [500, 500, 900, 300])
    
    strName4 = sprintf('Respiration #%d, PSD plots',str2num(data_num));
    saveas(h,strName4,'png');
    
    %% get cardiac signal slightly longer than BOLD so that spectrogram times are aligned to fMRI times
    ws = 2^14; % spectrogram window length. Needs to be specified beforehand because of specific matching procedure.
    fs = samples_per_second; %blood sampling rate
    fs2 = samples_per_second2; %respiration sampling rate
    
    clear beg_long fin_long beg_long2 fin_long2;
    
    if beg(chan_card) > ws/2 - TR*fs/2
        beg_long = beg(chan_card) - ws/2 + TR*fs/2;
    else
        beg_long = 1;
        %         disp(['Run ',num2str(run),' - signal too short at the beginning. Alignment wont be perfect.'])
    end
    if fin(chan_card) + ws/2 - TR*fs/2 < length(biopac_data.channels{chan_card}.data)
        fin_long = fin(chan_card) + ws/2 - TR*fs/2;
    else
        fin_long = length(biopac_data.channels{chan_card}.data);
        %         disp(['Run ',num2str(run),' - signal too short at the end. Alignment wont be perfect.'])
    end
    cardiac_long{exp_case_num} = biopac_data.channels{chan_card}.data(beg_long:fin_long)';
    
    %removing mean values
    cardiac_match{exp_case_num} = cardiac_match{exp_case_num} - mean(cardiac_match{exp_case_num},1);
    cardiac_long{exp_case_num} = cardiac_long{exp_case_num} - mean(cardiac_long{exp_case_num},1);
    
    %%respiration
    if beg(chan_resp) > ws/2 - TR*fs2/2
        beg_long2 = beg(chan_resp) - ws/2 + TR*fs2/2;
    else
        beg_long2 = 1;
        %         disp(['Run ',num2str(run),' - signal too short at the beginning. Alignment wont be perfect.'])
    end
    if fin(chan_resp) + ws/2 - TR*fs2/2 < length(biopac_data.channels{chan_resp}.data)
        fin_long2 = fin(chan_resp) + ws/2 - TR*fs2/2;
    else
        fin_long2 = length(biopac_data.channels{chan_resp}.data);
        %         disp(['Run ',num2str(run),' - signal too short at the end. Alignment wont be perfect.'])
    end
    respiratory_long{exp_case_num} = biopac_data.channels{chan_resp}.data(beg_long2:fin_long2)';
    
    %removing mean values
    respiratory_match{exp_case_num}     = respiratory_match{exp_case_num} - mean(respiratory_match{exp_case_num},1);
    respiratory_long{exp_case_num}      = respiratory_long{exp_case_num}  - mean(respiratory_long{exp_case_num},1);
% end
 
%% generate time
t_fmri = [0 : TR : (length(fmri_dummy)-1)*TR]';
t_cardiac = [0 : 1/fs : (length(cardiac_match{1})-1)/fs]';
t_respiratory = [0 : 1/fs2 : (length(respiratory_match{1})-1)/fs2]';

    
%% regress polynomial of degree n from fMRI data
if RegressPoly_On == 1
    for is=1:Total_Slice
%         for ir = 1 : size(Exp_Number_Array1,2)
            for iv = 1:size(cortex{is}{exp_case_num}, 2)
                p = polyfit(t_fmri, cortex{is}{exp_case_num}(:, iv), 3);
                reg_cortex{is}{exp_case_num}(:, iv) = cortex{is}{exp_case_num}(:, iv) - polyval(p, t_fmri);
            end
%         end
    end
    cortex = reg_cortex;
%     cortex{num_slice}{exp_case_num} = abs(squeeze(linescanning_data(:,:,num_slice)))';
            
    % just for checking the line-scanning signals
    for num_slice=1:Total_Slice
        h= figure();
        if num_slice == 1
            if num_slice == Total_Slice
                plot(1:size(linescanning_data,2),zscore(mean(abs(cortex{num_slice}{exp_case_num}'),1),0,2),'-k','LineWidth',plot_linewidth);
            else
                plot(1:size(linescanning_data,2),zscore(mean(abs(cortex{num_slice}{exp_case_num}'),1),0,2),'-b','LineWidth',plot_linewidth);
            end
        elseif num_slice == 2
            plot(1:size(linescanning_data,2),zscore(mean(abs(cortex{num_slice}{exp_case_num}'),1),0,2),'-r','LineWidth',plot_linewidth);
        else
            plot(1:size(linescanning_data,2),zscore(mean(abs(cortex{num_slice}{exp_case_num}'),1),0,2),'-k','LineWidth',plot_linewidth);
        end
        strName1 = sprintf('Slice #%d, RegressPolynomial averaged voxel time course',num_slice);
        xlim([0 640*(1/TR)])
        xticks(0:160*(1/TR):640*(1/TR))
        xticklabels({'0s','160s','320s','480s','640s'})
        %     xlabel('Time(sec)');
        ylim([-4 4])
        ylabel('fMRI Signal (a.u)');
        title(strName1);
        set(gcf,'color','w');
        set(gca,'FontSize',font_size,'FontWeight','Bold');
        set(gca,'linewidth',axs_line_width);
        set(gcf, 'Position', [500, 500, 900, 300])
        title(strName1);
        saveas(h,strName1,'png');
    end
end


%% regress Physiological Noise from fMRI data
%Following Hu, the cardiac phase is defined as Phi_t = 2*pi(t - t1)/(t2-t1)
%where t1 is the time of the R-wave peak in the cardiac cycle just preceding t, 
%and t2 is the time for the subsequent R-wave peak. 
if RegressPhysio_On == 1
   
    % blood_press_temp is original cardiac data
    %find the RR peak 
%     clear blood_press_temp;
%     blood_press_temp = synchorized_blood_press_temp';
    np = size(blood_press_temp,1);
%     blood_press_temp = blood_press_temp - mean(blood_press_temp);
    blood_press_temp = (blood_press_temp - mean(blood_press_temp))*(-1);%to find RR peak, see the paper, Cuff-Free Blood Pressure Estimation Using Pulse Transit Time and Heart Rate
    [temp_cardiac_envelop,temp_locs] = findpeaks(blood_press_temp(1:np,1),'MinPeakDistance',500);
    %downsampling carse interval is 2
%     [temp_cardiac_envelop,temp_locs] = findpeaks(blood_press_temp(1:np,1),'MinPeakDistance',2);
     h=figure();
    %     subplot(2,1,1);
    plot(1:size(blood_press_temp(1:np,1),1),blood_press_temp(1:np,1)*(-1),'r-');
    title('Original blood press')
    
    %     subplot(2,1,2);
    hold on
    plot(temp_locs,temp_cardiac_envelop(:,1)*(-1),'b*');
    title('Lower peak envelope blood pressure')
    ylim([-20 20])
    set(gcf,'color','w');
    set(gca,'FontSize',font_size,'FontWeight','Bold');
    set(gca,'linewidth',axs_line_width);
    set(gcf, 'Position', [500, 500, 1800, 500])
    if RS == 1
        strName3 = sprintf('Rest #%d, original biopac with peak envelope',str2num(data_num));
    else
        strName3 =  sprintf('Evoked #%d, original biopac with peak envelope',str2num(data_num));
    end
    legend('signal','peak envelope','Location','NorthEast','NumColumns',2)
    legend boxoff;
    xlim([0 np]);
    set(gcf,'color','w');
    saveas(h,strName3,'png');
     %% Good RR peak without Resampling
    cardiac_phase_time = 1:TR/(1/samples_per_second):np;
    nPhase = size(cardiac_phase_time,2);
%     nPeak = size(temp_locs,1);
    for phase_order = 1 : nPhase
        if cardiac_phase_time(phase_order) <= temp_locs(1) %copy just the preceding cycle
            one_R_peak = 1;
            next_R_peak = one_R_peak +1;
            one_R_peak_time  =  temp_locs(one_R_peak);
            next_R_peak_time  =  temp_locs(next_R_peak);
            zero_R_peak_time = one_R_peak_time - (next_R_peak_time - one_R_peak_time);
            %         for phase_order = 1 : one_R_peak_time
            Phi_card_time(phase_order) = 2*pi*(cardiac_phase_time(phase_order) - zero_R_peak_time)./(next_R_peak_time - one_R_peak_time);
        elseif cardiac_phase_time(phase_order) >= temp_locs(end) %copy just the preceding cycle
            one_R_peak_time = temp_locs(end-1);
            next_R_peak_time = temp_locs(end);
            Phi_card_time(phase_order) = 2*pi*(cardiac_phase_time(phase_order) - next_R_peak_time)./(next_R_peak_time - one_R_peak_time);
        else%     phase_order = 3;
            t_start = find(temp_locs < cardiac_phase_time(phase_order), 1, 'last');
            one_R_peak_time = temp_locs(t_start);
            next_R_peak_time = temp_locs(t_start+1);
            Phi_card_time(phase_order) = 2*pi*(cardiac_phase_time(phase_order) - one_R_peak_time)./(next_R_peak_time - one_R_peak_time);
        end
    end
   figure,
   plot(cos(Phi_card_time))
   %% Cardiac Fourier Expansion
   m_order = 2;
    for m = 1 : m_order
%         basis_cos_m_temp(:,m) = cos(m*Phi_card_time)';
%         basis_sin_m_temp(:,m) = sin(m*Phi_card_time)';
%         
%         basis_cos_m(:,m)  = resample(basis_cos_m_temp(:,m),Fs,samples_per_second);
%         basis_sin_m(:,m)  = resample(basis_sin_m_temp(:,m),Fs,samples_per_second);
        basis_cos_m(:,m) = cos(m*Phi_card_time)';
        basis_sin_m(:,m) = sin(m*Phi_card_time)';
        
%         basis_cos_m(:,m)  = resample(basis_cos_m_temp(:,m),Fs,samples_per_second);
%         basis_sin_m(:,m)  = resample(basis_sin_m_temp(:,m),Fs,samples_per_second);
    end
    
    for is=1:Total_Slice
       for iv = 1:size(cortex{is}{exp_case_num}, 2)
         
%                a_coeff(m) = regress(cortex{is}{exp_case_num}(:,iv), basis_cos_m(:,m));
%                b_coeff(m) = regress(cortex{is}{exp_case_num}(:,iv), basis_sin_m(:,m));
               
               a_coeff(:,iv) = regress(cortex{is}{exp_case_num}(:,iv), basis_cos_m(:,:));
               b_coeff(:,iv) = regress(cortex{is}{exp_case_num}(:,iv), basis_sin_m(:,:));
%                a_coeff(m) = regress(basis_cos_m(:,m),cortex{is}{exp_case_num}(:,iv));
%                b_coeff(m) = regress(basis_sin_m(:,m),cortex{is}{exp_case_num}(:,iv));
%                 a_coeff(m) = regress(blood_press_temp, basis_cos_m(:,m));
%                b_coeff(m) = regress(blood_press_temp, basis_sin_m(:,m));
%                a_coeff(m) = pinv(basis_cos_m(:,m))*cortex{is}{exp_case_num}(:, iv);
%                b_coeff(m) =
%                pinv(basis_sin_m(:,m))*cortex{is}{exp_case_num}(:, iv);
            for m = 1 : m_order
                cardiac_regress_match{exp_case_num}(:,m,iv) = a_coeff(m,iv)*basis_cos_m(:,m) + b_coeff(m,iv)*basis_sin_m(:,m);
%                 cardiac_regress_match{exp_case_num}(:,m,iv) = a_coeff(m,iv)* (basis_cos_m(:,m) + basis_sin_m(:,m));
%                a_coeff(m) = regress(cortex{is}{exp_case_num}(:, :), basis_cos_m(:,m));
%                b_coeff(m) = regress(cortex{is}{exp_case_num}(:, :), basis_sin_m(:,m));
%                cardiac_regress_match{exp_case_num}(:,m) = a_coeff(m)*basis_cos_m(:,m) + b_coeff(m)*basis_sin_m(:,m);
           end
%            for iv = 1:size(cortex{is}{exp_case_num}, 2)
%             w =5;
%            reg_card_cortex{is}{exp_case_num}(:, iv) = cortex{is}{exp_case_num}(:, iv) - squeeze(sum(cardiac_regress_match{exp_case_num}(:,:,iv),2));
           reg_card_cortex{is}{exp_case_num}(:, iv) = (cortex{is}{exp_case_num}(:, iv)) - squeeze(sum(cardiac_regress_match{exp_case_num}(:,:,iv),2));
%                 temp_card(:, iv, m) =  cardiac_regress_match{exp_case_num}(:,m,iv);
           card_only_cortex{is}{exp_case_num}(:, iv) =  squeeze(sum(cardiac_regress_match{exp_case_num}(:,:,iv),2));
%            end
       end
 %%
       h=figure();
       imagesc(card_only_cortex{is}{exp_case_num}')
       colormap(jet)
       colorbar;
        Xlim = get(gca, 'xlim');
       set(gca, 'XTick', linspace(Xlim(1), Xlim(end), 5));
       if TA == 640
           xticklabels({'0s','160s','320s','480s','640s'})
       elseif TA == 768
           xticklabels({'0s','192s','384s','576s','768s'})
       end
       ylabel('# of voxels');
       set(gcf,'color','w');
       set(gca,'FontSize',font_size,'FontWeight','Bold');
       set(gca,'linewidth',axs_line_width);
       strName20 = 'Cardiac spatial pattern';
       title(strName20);
       saveas(h,strName20,'png');
       %%
       res_card(:,:) = squeeze(mean(abs(card_only_cortex{is}{exp_case_num}'),1));
       %         res_line2 = zscore(abs(squeeze(res_line(:,:,num_slice))),0,2);
       res_card2 = res_card(:,:) - mean(res_card(:,:),2);
       %
       h=figure();
       
       %         [s, f, t, power] = spectrogram(res_line2,hanning(128),96,1024,Fs,'yaxis','power');
       [s_card, f_card, t_card, card_power(:,:)] = spectrogram(res_card2,hanning(128),96,1024,Fs,'yaxis','power');
       
       add_offset = min(min(10*log10(abs(card_power(:,:)))))*(-1)/2;
       surf(t_card, f_card, 10*log10(abs(card_power(:,:)))+add_offset, 'EdgeColor', 'none');
       axis xy;
       axis tight;
%        if Filter_On == 0
           colorbar_max = 30;
           colorbar_min = -20;
%        else
%            colorbar_max = 100;
%            colorbar_min = -50;
%        end
       colormap(jet);
       view(0,90);
       %     xlabel('Time(sec)');
       Xlim = get(gca, 'xlim');
       set(gca, 'XTick', linspace(Xlim(1), Xlim(end), 5));
       if TA == 640
           xticklabels({'0s','160s','320s','480s','640s'})
       elseif TA == 768
           xticklabels({'0s','192s','384s','576s','768s'})
       end
       colorbar;
       ylabel('Frequency(Hz)');
       set(gcf,'color','w');
       set(gca,'FontSize',font_size,'FontWeight','Bold');
       set(gca,'linewidth',axs_line_width);
       colormap(jet);
       c =colorbar('YTickLabel',{num2str(colorbar_min), num2str(colorbar_max)},'YTick', [colorbar_min colorbar_max]);
       c.Label.String = 'Power (dB)';
       caxis([colorbar_min colorbar_max])
       set(gcf, 'Position', [500, 500, 900, 500])
       %     strName3 = 'Blood Pressure Spectrogram (dB)';
       set(c,'YTickMode','manual');
       strName4 = sprintf('RETRICOR Cardiac Regressor Spectrogram');
       title(strName4);
       strName40 = sprintf('Line-scanning#%d,RETROICOR Cardiac Regressor Spectrogram(dB)',str2num(data_num));
       %         strName40 = strcat('Line-scanning', '#', strdata_num,',',strName4);
       saveas(h,strName40,'png');
       
       %% regressed cardiac
%        res_line(:,:,is) = squeeze(mean(abs(reg_card_cortex{is}{exp_case_num}'),1));
       res_line(:,:,is) = squeeze(mean(abs(reg_card_cortex{is}{exp_case_num}'),1));
       %         res_line2 = zscore(abs(squeeze(res_line(:,:,num_slice))),0,2);
       res_line2 = res_line(:,:,is) - mean(res_line(:,:,is),2);
       
       h=figure();
       [s_line, f_line, t_line, line_power(:,:,is)] = spectrogram(res_line2,hanning(128),96,1024,Fs,'yaxis','power');
       
       add_offset = min(min(10*log10(abs(line_power(:,:,is)))))*(-1)/2;
       surf(t_line, f_line, 10*log10(abs(line_power(:,:,is)))+add_offset, 'EdgeColor', 'none');
       axis xy;
       axis tight;
%        if Filter_On == 0
           colorbar_max = 50;
           colorbar_min = 0;
%        else
%            colorbar_max = 100;
%            colorbar_min = -50;
%        end
       colormap(jet);
       view(0,90);
       %     xlabel('Time(sec)');
       Xlim = get(gca, 'xlim');
       set(gca, 'XTick', linspace(Xlim(1), Xlim(end), 5));
       if TA == 640
           xticklabels({'0s','160s','320s','480s','640s'})
       elseif TA == 768
           xticklabels({'0s','192s','384s','576s','768s'})
       end
       colorbar;
       ylabel('Frequency(Hz)');
       set(gcf,'color','w');
       set(gca,'FontSize',font_size,'FontWeight','Bold');
       set(gca,'linewidth',axs_line_width);
       colormap(jet);
       c =colorbar('YTickLabel',{num2str(colorbar_min), num2str(colorbar_max)},'YTick', [colorbar_min colorbar_max]);
       c.Label.String = 'Power (dB)';
       caxis([colorbar_min colorbar_max])
       set(gcf, 'Position', [500, 500, 900, 500])
       %     strName3 = 'Blood Pressure Spectrogram (dB)';
       set(c,'YTickMode','manual');
       strName4 = sprintf('Slice #%d, Cardiac Regressed Line-scanning fMRI Spectrogram',is);
       title(strName4);
       strName40 = sprintf('Line-scanning#%d,Slice #%d,Cardiac Regressed Line-scanning fMRI Spectrogram(dB)',str2num(data_num),is);
       %         strName40 = strcat('Line-scanning', '#', strdata_num,',',strName4);
       saveas(h,strName40,'png');
       
       %% pure cardiac signal
%        res_line(:,:,is) = squeeze(mean(abs(reg_card_cortex{is}{exp_case_num}'),1));
%        res_cardiac(:,:) = squeeze((abs(cardiac_match{exp_case_num}')));
%        %         res_line2 = zscore(abs(squeeze(res_line(:,:,num_slice))),0,2);
%        res_cardiac2 = res_cardiac(:,:)- mean(res_cardiac(:,:),2);
%        
%        h=figure();
%        [s_line, f_line, t_line, card_power(:,:)] = spectrogram(res_cardiac2,hanning(128),96,1024,Fs,'yaxis','power');
%        
%        add_offset = min(min(10*log10(abs(card_power(:,:)))))*(-1)/2;
%        surf(t_line, f_line, 10*log10(abs(card_power(:,:)))+add_offset, 'EdgeColor', 'none');
%        axis xy;
%        axis tight;
% %        if Filter_On == 0
%            colorbar_max = 50;
%            colorbar_min = 0;
% %        else
% %            colorbar_max = 100;
% %            colorbar_min = -50;
% %        end
%        colormap(jet);
%        view(0,90);
%        %     xlabel('Time(sec)');
%        Xlim = get(gca, 'xlim');
%        set(gca, 'XTick', linspace(Xlim(1), Xlim(end), 5));
%        if TA == 640
%            xticklabels({'0s','160s','320s','480s','640s'})
%        elseif TA == 768
%            xticklabels({'0s','192s','384s','576s','768s'})
%        end
%        colorbar;
%        ylabel('Frequency(Hz)');
%        set(gcf,'color','w');
%        set(gca,'FontSize',font_size,'FontWeight','Bold');
%        set(gca,'linewidth',axs_line_width);
%        colormap(jet);
%        c =colorbar('YTickLabel',{num2str(colorbar_min), num2str(colorbar_max)},'YTick', [colorbar_min colorbar_max]);
%        c.Label.String = 'Power (dB)';
%        caxis([colorbar_min colorbar_max])
%        set(gcf, 'Position', [500, 500, 900, 500])
%        %     strName3 = 'Blood Pressure Spectrogram (dB)';
%        set(c,'YTickMode','manual');
%        strName4 = sprintf('Slice #%d, Blood Pressure Spectrogram',is);
%        title(strName4);
%        strName40 = sprintf('Line-scanning#%d,Blood Pressure Spectrogram(dB)',str2num(data_num),is);
%        %         strName40 = strcat('Line-scanning', '#', strdata_num,',',strName4);
%        saveas(h,strName40,'png');
       %% Line-scanning fMRI Powe spectrogram
       res_line(:,:,is) = squeeze(mean(abs(cortex{is}{exp_case_num}'),1));
       %         res_line2 = zscore(abs(squeeze(res_line(:,:,num_slice))),0,2);
       res_line2 = res_line(:,:,is) - mean(res_line(:,:,is),2);
       
       h=figure();
       [s_line, f_line, t_line, line_power(:,:,is)] = spectrogram(res_line2,hanning(128),96,1024,Fs,'yaxis','power');
       
       add_offset = min(min(10*log10(abs(line_power(:,:,is)))))*(-1)/2;
       surf(t_line, f_line, 10*log10(abs(line_power(:,:,is)))+add_offset, 'EdgeColor', 'none');
       axis xy;
       axis tight;
%        if Filter_On == 0
           colorbar_max = 50;
           colorbar_min = 0;
%        else
%            colorbar_max = 100;
%            colorbar_min = -50;
%        end
       colormap(jet);
       view(0,90);
       %     xlabel('Time(sec)');
       Xlim = get(gca, 'xlim');
       set(gca, 'XTick', linspace(Xlim(1), Xlim(end), 5));
       if TA == 640
           xticklabels({'0s','160s','320s','480s','640s'})
       elseif TA == 768
           xticklabels({'0s','192s','384s','576s','768s'})
       end
       yticklabels('manual')
       Ylim = get(gca, 'ylim');
       if Fs == 10
           set(gca, 'YTick', linspace(Ylim(1), Ylim(end), 6));
           yticklabels({'0','1','2','3','4','5'})
       elseif Fs == 20
           set(gca, 'YTick', linspace(Ylim(1), Ylim(end),11));
           yticklabels({'0','1','2','3','4','5','6','7','8','9','10'})
       end
    
       colorbar;
       ax1 = gca;
       ax1.YLabel.Visible = 'off';
       ax1.XLabel.Visible = 'off';
       set(gca,'xtick',[]);
       ylabel('Frequency(Hz)');
       set(gcf,'color','w');
       set(gca,'FontSize',font_size,'FontWeight','Bold');
       set(gca,'linewidth',axs_line_width);
       colormap(jet);
       c =colorbar('YTickLabel',{num2str(colorbar_min), num2str(colorbar_max)},'YTick', [colorbar_min colorbar_max]);
       c.Label.String = 'Power (dB)';
       caxis([colorbar_min colorbar_max])
       set(gcf, 'Position', [500, 500, 900, 600])
       %     strName3 = 'Blood Pressure Spectrogram (dB)';
       set(c,'YTickMode','manual');
       strName4 = sprintf('Slice #%d,Line-scanning fMRI Spectrogram',is);
%        title(strName4);
       if RS == 0
           strName40 = sprintf('Evoked #%d, Only Line-scanning fMRI Spectrogram(dB),Slice #%d',str2num(data_num),is);
       else
           strName40 = sprintf('Rest #%d, Only Line-scanning fMRI Spectrogram(dB),Slice #%d',str2num(data_num),is);
       end
       %         strName40 = strcat('Line-scanning', '#', strdata_num,',',strName4);
       saveas(h,strName40,'png');
       %% Zoomed line-scanning spectrgram
       ymax_freq = 1;
        h=figure();
       [s_line, f_line, t_line, line_power(:,:,is)] = spectrogram(res_line2,hanning(128),96,1024,Fs,'yaxis','power');
       
       add_offset = min(min(10*log10(abs(line_power(:,:,is)))))*(-1)/2;
       surf(t_line, f_line, 10*log10(abs(line_power(:,:,is)))+add_offset, 'EdgeColor', 'none');
       axis xy;
       axis tight;
%        if Filter_On == 0
           colorbar_max = 50;
           colorbar_min = 0;
%        else
%            colorbar_max = 100;
%            colorbar_min = -50;
%        end
       colormap(jet);
       view(0,90);
       %     xlabel('Time(sec)');
       Xlim = get(gca, 'xlim');
       set(gca, 'XTick', linspace(Xlim(1), Xlim(end), 5));
       if TA == 640
           xticklabels({'0s','160s','320s','480s','640s'})
       elseif TA == 768
           xticklabels({'0s','192s','384s','576s','768s'})
       end
       yticklabels('manual')
       Ylim = get(gca, 'ylim');
       if Fs == 10
           set(gca, 'YTick', linspace(Ylim(1), Ylim(end), 6));
           yticklabels({'0','1','2','3','4','5'})
       elseif Fs == 20
           set(gca, 'YTick', linspace(Ylim(1), Ylim(end),11));
           yticklabels({'0','1','2','3','4','5','6','7','8','9','10'})
       end
       ylim([0 ymax_freq]);
       colorbar;
       ax1 = gca;
       ax1.YLabel.Visible = 'off';
       ax1.XLabel.Visible = 'off';
       set(gca,'xtick',[]);
       ylabel('Frequency(Hz)');
       set(gcf,'color','w');
       set(gca,'FontSize',font_size,'FontWeight','Bold');
       set(gca,'linewidth',axs_line_width);
       colormap(jet);
       c =colorbar('YTickLabel',{num2str(colorbar_min), num2str(colorbar_max)},'YTick', [colorbar_min colorbar_max]);
       c.Label.String = 'Power (dB)';
       caxis([colorbar_min colorbar_max])
       set(gcf, 'Position', [500, 500, 900, 300])
       %     strName3 = 'Blood Pressure Spectrogram (dB)';
       set(c,'YTickMode','manual');
       strName4 = sprintf('Slice #%d,Line-scanning fMRI Spectrogram',is);
%        title(strName4);
       if RS == 0
           strName40 = sprintf('Evoked #%d, Zoomed Line-scanning fMRI Spectrogram(dB),Slice #%d',str2num(data_num),is);
       else
           strName40 = sprintf('Rest #%d, Zoomed Line-scanning fMRI Spectrogram(dB),Slice #%d',str2num(data_num),is);
       end
       %         strName40 = strcat('Line-scanning', '#', strdata_num,',',strName4);
       saveas(h,strName40,'png');
      
       %%
       %         x = zscore(sum(abs(cortical_depth_map(:,:)),1));
       x =res_line2;
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
       set(gcf, 'Position', [400, 400, 900, 300])
       set(gca,'FontSize',font_size,'FontWeight','Bold');
       set(gca,'linewidth',axs_line_width);
       set(gcf,'color','w');
       strName2 = sprintf('Slice #%d, Cardiac Regressed Averaged PSD with whole time course',is);
       title(strName2);
       saveas(h,strName2,'png');
    end
    %%
    
% for is=1:Total_Slice
% %         for ir = 1 : size(Exp_Number_Array1,2)
%             for iv = 1:size(cortex{is}{exp_case_num}, 2)
%                 p = polyfit(t_fmri, cortex{is}{exp_case_num}(:, iv), 3);
%                 reg_cortex{is}{exp_case_num}(:, iv) = cortex{is}{exp_case_num}(:, iv) - polyval(p, t_fmri);
%             end
% %         end
%     end
%     

    %back to original variable
%     cortex = reg_card_cortex;
            
    % just for checking the line-scanning signals
    for num_slice=1:Total_Slice
        h= figure();
        if num_slice == 1
            if num_slice == Total_Slice
                plot(1:size(reg_card_cortex{num_slice}{exp_case_num},1),zscore(mean(abs(reg_card_cortex{num_slice}{exp_case_num}'),1),0,2),'-k','LineWidth',plot_linewidth);
            else
                plot(1:size(reg_card_cortex{num_slice}{exp_case_num},1),zscore(mean(abs(reg_card_cortex{num_slice}{exp_case_num}'),1),0,2),'-b','LineWidth',plot_linewidth);
            end
        elseif num_slice == 2
            plot(1:size(reg_card_cortex{num_slice}{exp_case_num},1),zscore(mean(abs(reg_card_cortex{num_slice}{exp_case_num}'),1),0,2),'-r','LineWidth',plot_linewidth);
        else
            plot(1:size(reg_card_cortex{num_slice}{exp_case_num},1),zscore(mean(abs(reg_card_cortex{num_slice}{exp_case_num}'),1),0,2),'-k','LineWidth',plot_linewidth);
        end
        strName1 = sprintf('Slice #%d, Regressed Cardiac averaged voxel time course',num_slice);
        xlim([0 640*(1/TR)])
        xticks(0:160*(1/TR):640*(1/TR))
        xticklabels({'0s','160s','320s','480s','640s'})
        %     xlabel('Time(sec)');
        ylim([-4 4])
        ylabel('fMRI Signal (a.u)');
        title(strName1);
        set(gcf,'color','w');
        set(gca,'FontSize',font_size,'FontWeight','Bold');
        set(gca,'linewidth',axs_line_width);
        set(gcf, 'Position', [500, 500, 900, 300])
        title(strName1);
        saveas(h,strName1,'png');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Respiration Regressor, RETROICOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        synchorized_respiratory_press_temp (downsampled
%        interval = 15;
%        respiratory_press_temp = synchorized_respiratory_press_temp' - mean(synchorized_respiratory_press_temp');
       respiratory_press_temp = repiratory_press_temp-mean(repiratory_press_temp);
       np2 = size(respiratory_press_temp,1);
%        [temp_resp_envelop,~] = envelope(respiratory_press_temp(1,1:500), np2,'peak');
       %normalize 
       min_resp_val = min(respiratory_press_temp);
%        max_resp_val = max(respiratory_press_temp);
%        normal_respiratory_press = (respiratory_press_temp-min_resp_val)./(max_resp_val-min_resp_val);
       normal_respiratory_press = (respiratory_press_temp-min_resp_val);
       max_resp_val = max(normal_respiratory_press);
       min_resp_val2 = min(normal_respiratory_press); % equal to 0
       figure,
       plot(normal_respiratory_press(1:500,:),'-b');
       set(gcf,'color','w');
       set(gcf, 'Position', [500, 500, 900, 300])
       figure,
       plot(synchorized_respiratory_press_temp(1,1:500),'-b');
       set(gcf,'color','w');
       set(gcf, 'Position', [500, 500, 900, 300])
       
       h=figure();
       nbins = 100;
       histo = histogram(normal_respiratory_press,nbins);
       histo.FaceColor = [0.1 0.2 0.9];
       histo_respiratory_press = histo.Values';
       ylabel('Frequency');
       xlabel('Normalized respiration amplitude');
       set(gcf,'color','w');
       set(gca,'FontSize',15,'FontWeight','Bold');
       set(gca,'linewidth',1.5);
       set(gcf, 'Position', [500, 500, 900, 300])
       title('Histogram of Normalized Respiration');
       strName40 = sprintf('Line-scanning #%d, Histogram Respiration',str2num(data_num));
       saveas(h,strName40,'png');
       
       figure,
       plot(cumsum(histo_respiratory_press)/sum(histo_respiratory_press))
       %%
       %equalization transfer function
%        [J,T] = histeq(normal_respiratory_press);
%        figure
%        plot((0:255)/255,T);
%         sum = 0;
%         for idx = 0 : nbins
%             sum = sum + histo_respiratory_press(idx);
%             sum_resp_Histogram(idx) = sum;
%         end
       %%
       ntemp = np2;
       interval = 1000;
%        interval = 10;
       [temp_high_resp_envelop,temp_high_resp_locs] = findpeaks(normal_respiratory_press(1:ntemp,1),'MinPeakDistance',interval);
       [temp_low_resp_envelop,temp_low_resp_locs] = findpeaks(normal_respiratory_press(1:ntemp,1)*(-1),'MinPeakDistance',interval);
        h=figure();
    %     subplot(2,1,1);
        plot(1:size(normal_respiratory_press(1:ntemp,1),1),normal_respiratory_press(1:ntemp,1),'b-');

        %     subplot(2,1,2);
        hold on
        plot(temp_high_resp_locs,temp_high_resp_envelop(:,1),'k*');
        
%         hold on
%         plot(temp_low_resp_locs,temp_low_resp_envelop(:,1)*(-1),'g*');
        hold off
        title('Higher peak envelope respiratory freq')
        ylim([-5 20])
        xlim([0 np2])
        set(gcf,'color','w');
        set(gca,'FontSize',font_size,'FontWeight','Bold');
        set(gca,'linewidth',axs_line_width);
        set(gcf, 'Position', [500, 500, 1800, 500])
        if RS == 1
            strName3 = sprintf('Rest #%d, original biopac with respiration',str2num(data_num));
        else
            strName3 =  sprintf('Evoked #%d, original biopac with respiration',str2num(data_num));
        end
        legend('signal','peak envelope','Location','NorthEast','NumColumns',2)
        legend boxoff;
        set(gcf,'color','w');
        saveas(h,strName3,'png');
       %%
       respiratory_phase_time = round(1:(samples_per_second2*TR):np2);
%        respiratory_phase_time(1) = 1;
%        respiratory_phase_time = 1:np2;
       nPhase = size(respiratory_phase_time,2);
       
       for phase_order = 1 : nPhase
           % consider initial phase for sign +,-
           current_time = respiratory_phase_time(phase_order);
           normal_Amplitude_time = normal_respiratory_press(respiratory_phase_time(phase_order));
%            normal_Amplitude_time  = Amplitude_time-min_resp_val;
           
%            if isempty(find(current_time == temp_low_resp_locs)) ~= 1 %% if current exist on local min,max point
           if isempty(find(normal_Amplitude_time == min_resp_val2)) ~= 1 %% if current exist on local min,max point
               Phi_resp_time(phase_order) = 0;
           elseif isempty(find(normal_Amplitude_time == max_resp_val)) ~= 1 %% if current exist on local min,max point
               Phi_resp_time(phase_order) = pi;
           else
               if temp_high_resp_locs(1) < temp_low_resp_locs(1) && current_time < temp_high_resp_locs(1)
                   sign_dR_dt = 1;
               elseif temp_high_resp_locs(1) < temp_low_resp_locs(1) && current_time < temp_low_resp_locs(1) && current_time > temp_high_resp_locs(1)
                   sign_dR_dt = -1;
               elseif temp_high_resp_locs(1) > temp_low_resp_locs(1) && current_time < temp_low_resp_locs(1)
                   sign_dR_dt = -1;
               elseif temp_high_resp_locs(1) > temp_low_resp_locs(1) && current_time > temp_low_resp_locs(1) && current_time < temp_high_resp_locs(1) 
                   sign_dR_dt = 1;
               % consider end phase for sign +,-
               elseif temp_high_resp_locs(end) < temp_low_resp_locs(end) && current_time > temp_low_resp_locs(end)
                       sign_dR_dt = 1;
               elseif temp_high_resp_locs(end) < temp_low_resp_locs(end) && temp_high_resp_locs(end) < current_time && current_time < temp_low_resp_locs(end)
                       sign_dR_dt = -1;
               elseif temp_high_resp_locs(end) > temp_low_resp_locs(end) && current_time > temp_high_resp_locs(end)
                       sign_dR_dt = -1;
               elseif temp_high_resp_locs(end) > temp_low_resp_locs(end) && temp_low_resp_locs(end) < current_time && current_time < temp_high_resp_locs(end)
                       sign_dR_dt = 1;
               else
                   t_start = find(temp_high_resp_locs < current_time, 1, 'last');
                   one_Rt_high_time = temp_high_resp_locs(t_start);
                   next_Rt_high_time = temp_high_resp_locs(t_start+1);

                   t_lower = find(temp_low_resp_locs > one_Rt_high_time, 1, 'first');
                   one_Rt_low_time = temp_low_resp_locs(t_lower);
                   
                   if current_time > one_Rt_high_time && current_time < one_Rt_low_time %descending
                       sign_dR_dt = -1;
                   elseif current_time > one_Rt_low_time && current_time < next_Rt_high_time %ascending
                       sign_dR_dt = 1;
                   end
               end
               temp_sign(phase_order) = sign_dR_dt;
               Resp_bin_order(phase_order)  = ceil(normal_Amplitude_time/max_resp_val*100);
%                Phi_resp_time(phase_order) = sign_dR_dt * pi * (sum(histo_respiratory_press(1:Resp_bin_order(phase_order)))/sum(histo_respiratory_press));
               cum_sum_Resp = cumsum(histo_respiratory_press)/sum(histo_respiratory_press);
               pi_fraction = cum_sum_Resp(Resp_bin_order(phase_order));
               Phi_resp_time(phase_order) = sign_dR_dt * pi * pi_fraction;
           end
           
          
%            Phi_resp_time(phase_order) = sign_dR_dt * pi * (sum(histo_respiratory_press(1:Resp_bin_order(phase_order)))/sum(histo_respiratory_press));
%            Phi_resp_time(phase_order) =  pi * (sum(histo_respiratory_press(1:Resp_bin_order(phase_order)))/sum(histo_respiratory_press));
       end
       figure,plot(cum_sum_Resp);
%     for phase_order = 1 : nPhase
%         if find(sum(histo_respiratory_press(1:Resp_bin_order(phase_order))) ==0)
%             phase_order
%         end
%     end
    
    figure,plot(temp_sign,'*b')
    figure,plot(cos(Phi_resp_time))
    figure,plot(sin(Phi_resp_time))

    %%
    h=figure();
    nbins = 100;
    histo = histogram( normal_respiratory_press(respiratory_phase_time),nbins);
    histo.FaceColor = [0.1 0.2 0.9];
    histo_respiratory_press = histo.Values';
    ylabel('Frequency');
    xlabel('Normalized respiration amplitude');
    set(gcf,'color','w');
    set(gca,'FontSize',15,'FontWeight','Bold');
    set(gca,'linewidth',1.5);
    set(gcf, 'Position', [500, 500, 900, 300])
    title('Histogram of TR-sampled Normalized Respiration');
    strName40 = sprintf('Line-scanning #%d, TR-sampled Normalize Histogram Respiration',str2num(data_num));
    saveas(h,strName40,'png');
    
%%
        h=figure();
        plot(Phi_resp_time,'-k','LineWidth',1.5);
        hold on     
%         plot(respiratory_press_temp,'-b','LineWidth',1.2);
        plot(1:numel(respiratory_phase_time),respiratory_press_temp(respiratory_phase_time,1),'-b','LineWidth',1.5);
         hold off
        set(gcf,'color','w');
        set(gcf, 'Position', [500, 500, 900, 300])
        set(gca,'FontSize',font_size,'FontWeight','Bold');
        set(gca,'linewidth',axs_line_width);
        xlim([0 TA*(1/TR)])
        xticks(0:TA/4*(1/TR):TA*(1/TR))
        if TA == 640
            xticklabels({'0s','160s','320s','480s','640s'})
        elseif TA == 768
            xticklabels({'0s','192s','384s','576s','768s'})
        end
%         ylim([-pi pi])
        ylim([-pi-7 pi+7])
        xlim([1 100])
        yticks(-pi:pi:pi);
        yticklabels({'-π','0','π'});
        strName3 = 'Estimated Respiratory Phase';
        title(strName3);
        legend('Respiratory Phase', 'Respiratory Amp','Location','NorthEast','NumColumns',2)
        legend boxoff;
        strName31 = sprintf('Line-scanning#%d, Estimated Respiratory Phase',str2num(data_num));
        saveas(h,strName31,'png');
        
        
        %% Respiration Fourier Expansion
        clear basis2_cos_m
        clear basis2_sin_m
    m_order = 1;
    for m = 1 : m_order
%         basis_cos_m_temp(:,m) = cos(m*Phi_card_time)';
%         basis_sin_m_temp(:,m) = sin(m*Phi_card_time)';
%         
%         basis_cos_m(:,m)  = resample(basis_cos_m_temp(:,m),Fs,samples_per_second);
%         basis_sin_m(:,m)  = resample(basis_sin_m_temp(:,m),Fs,samples_per_second);
6
        basis2_cos_m(:,m) = cos(m*Phi_resp_time);
        basis2_sin_m(:,m) = sin(m*Phi_resp_time);
%          basis2_exp_m(:,m) = exp(m*Phi_resp_time);

%         basis2_cos_m_temp(:,m) = cos(m*Phi_resp_time)';
%         basis2_sin_m_temp(:,m) = sin(m*Phi_resp_time)';
%         
%         basis2_cos_m(:,m)  = resample(basis2_cos_m_temp(:,m),Fs,samples_per_second2);
%         basis2_sin_m(:,m)  = resample(basis2_sin_m_temp(:,m),Fs,samples_per_second2);
    end
   
    for is=1:Total_Slice
       for iv = 1:size(cortex{is}{exp_case_num}, 2)
           
%                a2_coeff(m) = regress(basis2_cos_m(:,m),cortex{is}{exp_case_num}(:,iv));
%                b2_coeff(m) = regress(basis2_sin_m(:,m),cortex{is}{exp_case_num}(:,iv));
               a2_coeff(:,iv) = regress(cortex{is}{exp_case_num}(:,iv), basis2_cos_m(:,:));
               b2_coeff(:,iv) = regress(cortex{is}{exp_case_num}(:,iv), basis2_sin_m(:,:));
%                a2_coeff(:,iv) = regress((cortex{is}{exp_case_num}(:,iv)), basis2_exp_m(:,:));

%                a2_coeff(m,iv) = pinv(basis2_cos_m(:,m))*cortex{is}{exp_case_num}(:,iv);
%                b2_coeff(m,iv) = pinv(basis2_sin_m(:,m))*cortex{is}{exp_case_num}(:,iv);
            for m = 1 : m_order
               respiration_regress_match{exp_case_num}(:,m,iv) = a2_coeff(m,iv)*basis2_cos_m(:,m) + b2_coeff(m,iv)*basis2_sin_m(:,m);
%                respiration_regress_match{exp_case_num}(:,m,iv) = a2_coeff(m,iv)* (basis2_exp_m(:,m));
            end
           reg_resp_cortex{is}{exp_case_num}(:, iv) = (cortex{is}{exp_case_num}(:, iv)) - squeeze(sum(respiration_regress_match{exp_case_num}(:,:,iv),2));
%            reg_resp_cortex{is}{exp_case_num}(:, iv) = (cortex{is}{exp_case_num}(:, iv));
%                 temp_card(:, iv, m) =  cardiac_regress_match{exp_case_num}(:,m,iv);
           resp_only_cortex{is}{exp_case_num}(:, iv) =  squeeze(sum(respiration_regress_match{exp_case_num}(:,:,iv),2));
       end
       figure,plot(a2_coeff)
%        figure,plot(b2_coeff)
 %%
       h=figure();
       imagesc(resp_only_cortex{is}{exp_case_num}')
       colormap(jet)
       colorbar;
        Xlim = get(gca, 'xlim');
       set(gca, 'XTick', linspace(Xlim(1), Xlim(end), 5));
       if TA == 640
           xticklabels({'0s','160s','320s','480s','640s'})
       elseif TA == 768
           xticklabels({'0s','192s','384s','576s','768s'})
       end
       ylabel('# of voxels');
       set(gcf,'color','w');
       set(gca,'FontSize',font_size,'FontWeight','Bold');
       set(gca,'linewidth',axs_line_width);
       strName20 = 'Respiration spatial pattern';
       title(strName20);
       saveas(h,strName20,'png');
       
       res_resp(:,:) = squeeze(mean(abs(resp_only_cortex{is}{exp_case_num}'),1));
       %         res_line2 = zscore(abs(squeeze(res_line(:,:,num_slice))),0,2);
       res_resp2 = res_resp(:,:) - mean(res_resp(:,:),2);
       
       h=figure();
       
       %         [s, f, t, power] = spectrogram(res_line2,hanning(128),96,1024,Fs,'yaxis','power');
       [s_resp, f_resp, t_resp, resp_power(:,:)] = spectrogram(res_resp2,hanning(128),96,1024,Fs,'yaxis','power');
       
       add_offset = min(min(10*log10(abs(resp_power(:,:)))))*(-1)/2;
       surf(t_resp, f_resp, 10*log10(abs(resp_power(:,:)))+add_offset, 'EdgeColor', 'none');
       axis xy;
       axis tight;
%        if Filter_On == 0
%            colorbar_max = 30;
%            colorbar_min = -20;
%        else
           colorbar_max = 50;
           colorbar_min = 0;
%        end
       colormap(jet);
       view(0,90);
       %     xlabel('Time(sec)');
       Xlim = get(gca, 'xlim');
       set(gca, 'XTick', linspace(Xlim(1), Xlim(end), 5));
       if TA == 640
           xticklabels({'0s','160s','320s','480s','640s'})
       elseif TA == 768
           xticklabels({'0s','192s','384s','576s','768s'})
       end
       colorbar;
       ylabel('Frequency(Hz)');
       set(gcf,'color','w');
       set(gca,'FontSize',font_size,'FontWeight','Bold');
       set(gca,'linewidth',axs_line_width);
       colormap(jet);
       c =colorbar('YTickLabel',{num2str(colorbar_min), num2str(colorbar_max)},'YTick', [colorbar_min colorbar_max]);
       c.Label.String = 'Power (dB)';
       caxis([colorbar_min colorbar_max])
       set(gcf, 'Position', [500, 500, 900, 500])
       %     strName3 = 'Blood Pressure Spectrogram (dB)';
       set(c,'YTickMode','manual');
       strName4 = sprintf('RETROICOR Respiration Regressor Spectrogram');
       title(strName4);
       strName40 = sprintf('Line-scanning#%d,RETRICOR Respiration Regressor Spectrogram(dB)',str2num(data_num));
       %         strName40 = strcat('Line-scanning', '#', strdata_num,',',strName4);
       saveas(h,strName40,'png');
       
       
       res_line(:,:,is) = squeeze(mean(abs(reg_resp_cortex{is}{exp_case_num}'),1));
       %         res_line2 = zscore(abs(squeeze(res_line(:,:,num_slice))),0,2);
       res_line2 = res_line(:,:,is) - mean(res_line(:,:,is),2);
       
       h=figure();
       [s_line, f_line, t_line, line_power(:,:,is)] = spectrogram(res_line2,hanning(128),96,1024,Fs,'yaxis','power');
       
       add_offset = min(min(10*log10(abs(line_power(:,:,is)))))*(-1)/2;
       surf(t_line, f_line, 10*log10(abs(line_power(:,:,is)))+add_offset, 'EdgeColor', 'none');
       axis xy;
       axis tight;
%        if Filter_On == 0
           colorbar_max = 50;
           colorbar_min = 0;
%        else
%            colorbar_max = 100;
%            colorbar_min = -50;
%        end
       colormap(jet);
       view(0,90);
       %     xlabel('Time(sec)');
       Xlim = get(gca, 'xlim');
       set(gca, 'XTick', linspace(Xlim(1), Xlim(end), 5));
       if TA == 640
           xticklabels({'0s','160s','320s','480s','640s'})
       elseif TA == 768
           xticklabels({'0s','192s','384s','576s','768s'})
       end
       colorbar;
       ylabel('Frequency(Hz)');
       set(gcf,'color','w');
       set(gca,'FontSize',font_size,'FontWeight','Bold');
       set(gca,'linewidth',axs_line_width);
       colormap(jet);
       c =colorbar('YTickLabel',{num2str(colorbar_min), num2str(colorbar_max)},'YTick', [colorbar_min colorbar_max]);
       c.Label.String = 'Power (dB)';
       caxis([colorbar_min colorbar_max])
       set(gcf, 'Position', [500, 500, 900, 500])
       %     strName3 = 'Blood Pressure Spectrogram (dB)';
       set(c,'YTickMode','manual');
       strName4 = sprintf('Slice #%d, Respiration Regressed Line-scanning fMRI Spectrogram',is);
       title(strName4);
       strName40 = sprintf('Line-scanning#%d,Slice #%d,Respiration Regressed Line-scanning fMRI Spectrogram(dB)',str2num(data_num),is);
       %         strName40 = strcat('Line-scanning', '#', strdata_num,',',strName4);
       saveas(h,strName40,'png');
       
       % just for checking the line-scanning signals
        h= figure();
        if is == 1
            if is == Total_Slice
                plot(1:size(reg_resp_cortex{is}{exp_case_num},1),zscore(mean(abs(reg_resp_cortex{is}{exp_case_num}'),1),0,2),'-k','LineWidth',plot_linewidth);
            else
                plot(1:size(reg_resp_cortex{is}{exp_case_num},1),zscore(mean(abs(reg_resp_cortex{is}{exp_case_num}'),1),0,2),'-b','LineWidth',plot_linewidth);
            end
        elseif is == 2
            plot(1:size(reg_resp_cortex{is}{exp_case_num},1),zscore(mean(abs(reg_resp_cortex{is}{exp_case_num}'),1),0,2),'-r','LineWidth',plot_linewidth);
        else
            plot(1:size(reg_resp_cortex{is}{exp_case_num},1),zscore(mean(abs(reg_resp_cortex{is}{exp_case_num}'),1),0,2),'-k','LineWidth',plot_linewidth);
        end
        strName1 = sprintf('Slice #%d, Regressed Respiration averaged voxel time course',is);
        xlim([0 640*(1/TR)])
        xticks(0:160*(1/TR):640*(1/TR))
        xticklabels({'0s','160s','320s','480s','640s'})
        %     xlabel('Time(sec)');
        ylim([-4 4])
        ylabel('fMRI Signal (a.u)');
        title(strName1);
        set(gcf,'color','w');
        set(gca,'FontSize',font_size,'FontWeight','Bold');
        set(gca,'linewidth',axs_line_width);
        set(gcf, 'Position', [500, 500, 900, 300])
        title(strName1);
        saveas(h,strName1,'png');
    
    end
    if Cardiac_Reg_On == 1 && Respiration_Reg_On == 0
        cortex = reg_card_cortex;
    elseif Cardiac_Reg_On == 0 && Respiration_Reg_On == 1
        cortex = reg_resp_cortex;
    elseif Cardiac_Reg_On == 1 && Respiration_Reg_On == 1
         for is=1:Total_Slice
             cortex{is}{exp_case_num}(:,:) = cortex{is}{exp_case_num}(:,:) - card_only_cortex{is}{exp_case_num}(:,:) - resp_only_cortex{is}{exp_case_num}(:,:);
         end
    end
end
%% fMRI - between slice lags
% for ir = 1:length(Exp_Number_Array1)
if Total_Slice == 3
    [xxx12(exp_case_num, :), lag] = xcorr(zscore(mean(cortex{1}{exp_case_num}, 2)), zscore(mean(cortex{2}{exp_case_num}, 2)), 100, 'coeff');
    [xxx13(exp_case_num, :), lag] = xcorr(zscore(mean(cortex{1}{exp_case_num}, 2)), zscore(mean(cortex{3}{exp_case_num}, 2)), 100, 'coeff');
    [xxx23(exp_case_num, :), lag] = xcorr(zscore(mean(cortex{2}{exp_case_num}, 2)), zscore(mean(cortex{3}{exp_case_num}, 2)), 100, 'coeff');
    h = figure; set(h, 'WindowStyle', 'docked', 'NumberTitle', 'off', 'Name', ['fmri xcorr s12'])
    plot(lag*TR, mean(xxx12, 1))
    h = figure; set(h, 'WindowStyle', 'docked', 'NumberTitle', 'off', 'Name', ['fmri xcorr s13'])
    plot(lag*TR, mean(xxx13, 1))
    h = figure; set(h, 'WindowStyle', 'docked', 'NumberTitle', 'off', 'Name', ['fmri xcorr s23'])
    plot(lag*TR, mean(xxx23, 1))
end
%% compute spectrograms & extract freq. band signals
disp('compute spectrograms for cardiac beat')
% for ir = 1:length(Exp_Number_Array1)
    clear S_cardiac F_cardiac T_cardiac Cardiac_Power;
   
    win = hanning(ws);
    ovrl = ws - floor((length(cardiac_long{exp_case_num})-ws)/(length(fmri_dummy)-1));
    nfft = 2^15;
    [S_cardiac,F_cardiac,T_cardiac,Cardiac_Power] = spectrogram(cardiac_long{exp_case_num}, win, ovrl, nfft, fs, 'yaxis','power');
    
    if exp_case_num <= n_plots
        f_int = F_cardiac < 10 & F_cardiac>0;
        h = figure; set(h, 'WindowStyle', 'docked', 'NumberTitle', 'off', 'Name', ['cardiac spectrogram r', num2str(Exp_Number_Array1(exp_case_num))])
        %         imagesc(T, F(f_int), log10(abs(S(f_int, :)))), box off
        imagesc(T_cardiac, F_cardiac(f_int), abs(S_cardiac(f_int, :))), box off
        ax = gca;
        ax.YDir = 'normal';
        title(['spectrogram, trial ', data_num])
        box off
    end
    
    %%
     %calculate freq. band powers of blood pressure
%     f_low = 5; f_high = 8;
%     f_low = 6.5; f_high = 7.5; %for 3slice, 09062018
    f_low = card_f_low; f_high = card_f_high; %for 1 slice, 09252018 
%     cardiac_long_69{ir} = mean(abs(S_cardiac(F_cardiac<=f_high & F_cardiac>=f_low, :)))';
    t_start = find(F_cardiac<=f_low,1, 'last');
    t_end = find(F_cardiac>=f_high,1, 'first');
    
    cardiac_freq_67{exp_case_num}  = mean(10*log10(Cardiac_Power(t_start:t_end,:)),1);
    %     blood_power_wods_10 = mean(10*log10((blood_power_wods(:,:))),1);
%     cardiac_freq_69{ir}  = mean(abs(S_cardiac(F_cardiac<=f_high & F_cardiac>=f_low, :)),1);
    h=figure();
    plot(1:size(cardiac_freq_67{exp_case_num} ,2), cardiac_freq_67{exp_case_num}(1,:) ,'m-','LineWidth',plot_linewidth);
    xlim([0 TA*(1/TR)])
    xticks(0:TA/4*(1/TR):TA*(1/TR))
    if TA == 640
        xticklabels({'0s','160s','320s','480s','640s'})
    elseif TA == 768
        xticklabels({'0s','192s','384s','576s','768s'})
    end
    set(gcf,'color','w');
    set(gcf, 'Position', [500, 500, 1200, 300])
    set(gca,'FontSize',font_size,'FontWeight','Bold');
    set(gca,'linewidth',axs_line_width)
    strName1 = sprintf('Blood Pressure, %d-%dHz Band Power Time Course',f_low,f_high);
%     title(strName1);
    strName30 = strcat(strName1, '#', data_num);
    saveas(h,strName30,'png');
    %% Blood Pressure Gradient
    
    h=figure();
    plot(1:size(cardiac_freq_67{exp_case_num} ,2), gradient(cardiac_freq_67{exp_case_num}(1,:)) ,'b-','LineWidth',plot_linewidth);
    xlim([0 TA*(1/TR)])
    xticks(0:TA/4*(1/TR):TA*(1/TR))
    if TA == 640
        xticklabels({'0s','160s','320s','480s','640s'})
    elseif TA == 768
        xticklabels({'0s','192s','384s','576s','768s'})
    end
    set(gcf,'color','w');
    set(gcf, 'Position', [500, 500, 900, 300])
    set(gca,'FontSize',font_size,'FontWeight','Bold');
    set(gca,'linewidth',axs_line_width)
    strName1 = sprintf('Blood Pressure, %d-%dHz Band Power 1st Gradient Time Course',f_low,f_high);
    title(strName1);
    strName30 = sprintf('Blood Pressure, %d-%dHz Band Power 1st Gradient Time Course #%d',floor(f_low),floor(f_high),str2num(data_num));
%     strName30 = strcat(strName1, '#', data_num);
    saveas(h,strName30,'png');
    %% Blood Pressure 2nd Gradient
    
    h=figure();
    plot(1:size(cardiac_freq_67{exp_case_num} ,2), gradient(gradient(cardiac_freq_67{exp_case_num}(1,:))) ,'b-','LineWidth',plot_linewidth);
    xlim([0 TA*(1/TR)])
    xticks(0:TA/4*(1/TR):TA*(1/TR))
    if TA == 640
        xticklabels({'0s','160s','320s','480s','640s'})
    elseif TA == 768
        xticklabels({'0s','192s','384s','576s','768s'})
    end
    set(gcf,'color','w');
    set(gcf, 'Position', [500, 500, 900, 300])
    set(gca,'FontSize',font_size,'FontWeight','Bold');
    set(gca,'linewidth',axs_line_width)
    strName1 = sprintf('Blood Pressure, %d-%dHz Band Power 2nd Gradient Time Course',f_low,f_high);
    title(strName1);
    strName30 = sprintf('Blood Pressure, %d-%dHz Band Power 2nd Gradient Time Course #%d',floor(f_low),floor(f_high),str2num(data_num));
%     strName30 = strcat(strName1, '#', data_num);
    saveas(h,strName30,'png');
% end
%%
disp('compute spectrograms for respiratory freq')
% for ir = 1:length(Exp_Number_Array1)
    clear S_respiratory F_respiratory T_respiratory Respiratory_Power;
    
    data_num = num2str(Exp_Number_Array1(exp_case_num));
    win = hanning(ws);
    ovrl = ws - floor((length(respiratory_long{exp_case_num})-ws)/(length(fmri_dummy)-1));
    nfft = 2^15;
    [S_respiratory,F_respiratory,T_respiratory, Respiratory_Power] = spectrogram(respiratory_long{exp_case_num}, win, ovrl, nfft, fs2, 'yaxis','power');
    
    if exp_case_num <= n_plots
        f_int = F_respiratory < 10 & F_respiratory>0;
        h = figure; set(h, 'WindowStyle', 'docked', 'NumberTitle', 'off', 'Name', ['resp spectrogram r', num2str(Exp_Number_Array1(exp_case_num))])
        %         imagesc(T, F(f_int), log10(abs(S(f_int, :)))), box off
        imagesc(T_respiratory, F_respiratory(f_int), abs(S_respiratory(f_int, :))), box off
        ax = gca;
        ax.YDir = 'normal';
        title(['spectrogram, trial ', data_num])
        box off
    end
    
    %%
     %calculate freq. band powers of blood pressure
%     f_low = 0; f_high = 10;
    f_low = 0.5; f_high = 1.5;
%     cardiac_long_69{ir} = mean(abs(S_cardiac(F_cardiac<=f_high & F_cardiac>=f_low, :)))';
    t_start = find(F_respiratory<=f_low,1, 'last');
    t_end = find(F_respiratory>=f_high,1, 'first');
    
    respiratory_freq_01{exp_case_num}  = mean(10*log10(Respiratory_Power(t_start:t_end,:)),1);
    %     blood_power_wods_10 = mean(10*log10((blood_power_wods(:,:))),1);
%     cardiac_freq_69{ir}  = mean(abs(S_cardiac(F_cardiac<=f_high & F_cardiac>=f_low, :)),1);
    h=figure();
    plot(1:size(respiratory_freq_01{exp_case_num} ,2), respiratory_freq_01{exp_case_num}(1,:) ,'g-','LineWidth',plot_linewidth);
    xlim([0 TA*(1/TR)])
    xticks(0:TA/4*(1/TR):TA*(1/TR))
    if TA == 640
        xticklabels({'0s','160s','320s','480s','640s'})
    elseif TA == 768
        xticklabels({'0s','192s','384s','576s','768s'})
    end
    set(gcf,'color','w');
    set(gcf, 'Position', [500, 500, 900, 300])
    set(gca,'FontSize',font_size,'FontWeight','Bold');
    set(gca,'linewidth',axs_line_width)
    strName1 = sprintf('Respiratory Freq, %1.1f-%1.1fHz Band Power Time Course',f_low,f_high);
    title(strName1);
    strName3 = sprintf('Respiratory Freq, %d-%dHz Band Power Time Course #%d',floor(f_low),floor(f_high),str2num(data_num));
    saveas(h,strName3,'png');
    
    %% Respiratory Freq Gradient
    
    h=figure();
    plot(1:size(respiratory_freq_01{exp_case_num} ,2), gradient(respiratory_freq_01{exp_case_num}(1,:)) ,'g-','LineWidth',plot_linewidth);
    xlim([0 TA*(1/TR)])
    xticks(0:TA/4*(1/TR):TA*(1/TR))
    if TA == 640
        xticklabels({'0s','160s','320s','480s','640s'})
    elseif TA == 768
        xticklabels({'0s','192s','384s','576s','768s'})
    end
    set(gcf,'color','w');
    set(gcf, 'Position', [500, 500, 900, 300])
    set(gca,'FontSize',font_size,'FontWeight','Bold');
    set(gca,'linewidth',axs_line_width)
    strName1 = sprintf('Respiratory Freq, %1.1f-%1.1fHz Band Power 1st Gradient Time Course',f_low,f_high);
    title(strName1);
    strName30 = sprintf('Respiratory Freq, %1.1f-%1.1fHz Band Power 1st Gradient Time Course #%d',(f_low),(f_high),str2num(data_num));
%     strName30 = strcat(strName1, '#', data_num);
    saveas(h,strName30,'png');
% end
%% with matched cardiac signal
% for ir = 1:length(Exp_Number_Array1)
    
    data_num = num2str(Exp_Number_Array1(exp_case_num));
    synchorized_blood_press = resample(cardiac_match{exp_case_num},Fs1,fs)';
%      synchorized_blood_press = sum(cardiac_regress_match{exp_case_num},2);
    res_blood2 = synchorized_blood_press';
    
    nfft2 = 2^10;
    ws2 = 2^7;
    ovrl2 = ws2-ws2*1/4;
    
    h=figure();
    
    [s_blood, f_blood, t_blood, syn_blood_power] = spectrogram(res_blood2,hanning(128),96,1024,Fs1,'yaxis','power');
%     [s, f, t, power] = spectrogram(res_x2,hanning(128),ovrl,nfft,Fs1,'yaxis','power');
%     [s_blood, f_blood, t_blood, blood_power] = spectrogram(res_blood2,hanning(ws2),ovrl2,nfft2,Fs1,'yaxis','power');
%     [s_wods, f_wods, t_wods, blood_power_wods] = spectrogram(res_wo_ds,hanning(ws2),ovrl2,nfft2,Fs1,'yaxis','power');
%     size(blood_power_wods)
    
    add_offset2 = min(min(10*log10(abs(syn_blood_power))))*(-1)/2;
    surf(t_blood, f_blood, 10*log10(abs(syn_blood_power))+add_offset2, 'EdgeColor', 'none');

    axis xy;
    axis tight;
    colorbar_max = 80;
    colorbar_min = 30;
    colormap(jet); 
    view(0,90);
%     xlabel('Time(sec)');
    Xlim = get(gca, 'xlim');
    set(gca, 'XTick', linspace(Xlim(1), Xlim(end), 5));
    if TA == 640
        xticklabels({'0s','160s','320s','480s','640s'})
    elseif TA == 768
        xticklabels({'0s','192s','384s','576s','768s'})
    end
    colorbar;
    ylabel('Frequency(Hz)');
    set(gcf,'color','w');
    set(gca,'FontSize',font_size,'FontWeight','Bold');
    set(gca,'linewidth',axs_line_width);
    colormap(jet);
    c =colorbar('YTickLabel',{num2str(colorbar_min), num2str(colorbar_max)},'YTick', [colorbar_min colorbar_max]); 
    c.Label.String = 'Power (dB)';
    caxis([colorbar_min colorbar_max])
    set(gcf, 'Position', [500, 500, 800, 500])
%     strName3 = 'Blood Pressure Spectrogram (dB)';
    set(c,'YTickMode','manual'); 
    strName4 = 'Blood Pressure Spectrogram';
    title(strName4);

    strName3 = sprintf('Blood Pressure #%d, Spectrogram (dB)',str2num(data_num));
    saveas(h,strName3,'png');
%     save res_blood2.mat 
% end

%% with matched respiratory frequency
% for ir = 1:length(Exp_Number_Array1)
    
    synchorized_respiratory_freq = resample(respiratory_match{exp_case_num},Fs1,fs2)';
    res_resp2 = synchorized_respiratory_freq';
    
    nfft2 = 2^10;
    ws2 = 2^7;
    ovrl2 = ws2-ws2*1/4;
    
    h=figure();
    
    [s_resp, f_resp, t_resp, syn_resp_power] = spectrogram(res_resp2,hanning(128),96,1024,Fs1,'yaxis','power');
%     [s, f, t, power] = spectrogram(res_x2,hanning(128),ovrl,nfft,Fs1,'yaxis','power');
%     [s_blood, f_blood, t_blood, blood_power] = spectrogram(res_blood2,hanning(ws2),ovrl2,nfft2,Fs1,'yaxis','power');
%     [s_wods, f_wods, t_wods, blood_power_wods] = spectrogram(res_wo_ds,hanning(ws2),ovrl2,nfft2,Fs1,'yaxis','power');
%     size(blood_power_wods)
    
    add_offset2 = min(min(10*log10(abs(syn_resp_power))))*(-1)/2;
    surf(t_resp, f_resp, 10*log10(abs(syn_resp_power))+add_offset2, 'EdgeColor', 'none');

    axis xy;
    axis tight;
    colorbar_max = 70;
    colorbar_min = 0;
    colormap(jet); 
    view(0,90);
%     xlabel('Time(sec)');
    Xlim = get(gca, 'xlim');
    set(gca, 'XTick', linspace(Xlim(1), Xlim(end), 5));
    if TA == 640
        xticklabels({'0s','160s','320s','480s','640s'})
    elseif TA == 768
        xticklabels({'0s','192s','384s','576s','768s'})
    end
    colorbar;
    ylabel('Frequency(Hz)');
    set(gcf,'color','w');
    set(gca,'FontSize',font_size,'FontWeight','Bold');
    set(gca,'linewidth',axs_line_width);
    colormap(jet);
    c =colorbar('YTickLabel',{num2str(colorbar_min), num2str(colorbar_max)},'YTick', [colorbar_min colorbar_max]); 
    c.Label.String = 'Power (dB)';
    caxis([colorbar_min colorbar_max])
    set(gcf, 'Position', [500, 500, 900, 500])
%     strName3 = 'Blood Pressure Spectrogram (dB)';
    set(c,'YTickMode','manual'); 
    strName4 = 'Respiratory Freq Spectrogram';
    title(strName4);

    strName3 = sprintf('Respiratory Freq #%d, Spectrogram (dB)',str2num(data_num));
    saveas(h,strName3,'png');
% end

        %%
%         linescanning_data
        
        for num_slice = 1 : Total_Slice
%             linescanning_data(:,:,num_slice) = cortex{num_slice}{exp_case_num}';
            res_line(:,:,num_slice) = squeeze(mean(abs(cortex{num_slice}{exp_case_num}'),1));
            %         res_line2 = zscore(abs(squeeze(res_line(:,:,num_slice))),0,2);
            res_line2 = res_line(:,:,num_slice) - mean(res_line(:,:,num_slice),2);
            
            h=figure();
            
            %         [s, f, t, power] = spectrogram(res_line2,hanning(128),96,1024,Fs,'yaxis','power');
            [s_line, f_line, t_line, line_power(:,:,num_slice)] = spectrogram(res_line2,hanning(128),96,1024,Fs,'yaxis','power');
            
            add_offset = min(min(10*log10(abs(line_power(:,:,num_slice)))))*(-1)/2;
            surf(t_line, f_line, 10*log10(abs(line_power(:,:,num_slice)))+add_offset, 'EdgeColor', 'none');
            axis xy;
            axis tight;
            if Filter_On == 0
                colorbar_max = 50;
                colorbar_min = 0;
            else
                colorbar_max = 100;
                colorbar_min = -50;
            end
            colormap(jet);
            view(0,90);
            %     xlabel('Time(sec)');
            Xlim = get(gca, 'xlim');
            set(gca, 'XTick', linspace(Xlim(1), Xlim(end), 5));
            if TA == 640
                xticklabels({'0s','160s','320s','480s','640s'})
            elseif TA == 768
                xticklabels({'0s','192s','384s','576s','768s'})
            end
            colorbar;
%             ylabel('Frequency(Hz)');
            set(gcf,'color','w');
            set(gca,'FontSize',font_size,'FontWeight','Bold');
            set(gca,'linewidth',axs_line_width);
            colormap(jet);
            c =colorbar('YTickLabel',{num2str(colorbar_min), num2str(colorbar_max)},'YTick', [colorbar_min colorbar_max]);
            c.Label.String = 'Power (dB)';
            caxis([colorbar_min colorbar_max])
            set(gcf, 'Position', [500, 500, 900, 500])
            %     strName3 = 'Blood Pressure Spectrogram (dB)';
            set(c,'YTickMode','manual');
            strName4 = sprintf('Slice #%d, Averaged voxel Spectrogram',num_slice);
%             title(strName4);
            strName40 = sprintf('Line-scanning #%d, Slice #%d, Spectrogram (dB)',str2num(data_num),num_slice);
            %         strName40 = strcat('Line-scanning', '#', strdata_num,',',strName4);
            saveas(h,strName40,'png');
            %%
            %         x = zscore(sum(abs(cortical_depth_map(:,:)),1));
            x =res_line2;
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
            set(gcf, 'Position', [400, 400, 900, 300])
            set(gca,'FontSize',font_size,'FontWeight','Bold');
            set(gca,'linewidth',axs_line_width);
            set(gcf,'color','w');
            strName2 = sprintf('Slice #%d, Averaged PSD with whole time course',num_slice);
            title(strName2);
            saveas(h,strName2,'png');
            
            
            %%
            %line-scanning data
            nCortical_Voxel = size(cortex{num_slice}{exp_case_num},2);
            if Spatial_Resolution == 0.05
                Average_voxel_layers(1,:)  = mean(cortex{num_slice}{exp_case_num}(:,1:3),2);%L1,0~0.15 mm, averaged voxel
                Average_voxel_layers(2,:)  = mean(cortex{num_slice}{exp_case_num}(:,4:11),2);%L2/3,0.15~0.55 mm, averaged voxel
                Average_voxel_layers(3,:)  = mean(cortex{num_slice}{exp_case_num}(:,12:16),2);%L4,0.55~0.8 mm, averaged voxel
                Average_voxel_layers(4,:)  = mean(cortex{num_slice}{exp_case_num}(:,17:27),2);%L5,0.8~1.4 mm, averaged voxel
                Average_voxel_layers(5,:)  = mean(cortex{num_slice}{exp_case_num}(:,28:40),2);%L6,1.4~2.0 mm, averaged voxel
                if nCortical_Voxel > 40
                    Average_voxel_layers(6,:)  = mean(cortex{num_slice}{exp_case_num}(:,41:nCortical_Voxel),2);%WM,2.1~3.0 mm, averaged voxel
                end
                
            elseif Spatial_Resolution == 0.1
                Average_voxel_layers(1,:)  = mean(cortex{num_slice}{exp_case_num}(:,1:2),2);%L1,0~0.15 mm, averaged voxel
                Average_voxel_layers(2,:)  = mean(cortex{num_slice}{exp_case_num}(:,3:5),2);%L2/3,0.15~0.55 mm, averaged voxel
                Average_voxel_layers(3,:)  = mean(cortex{num_slice}{exp_case_num}(:,6:8),2);%L4,0.55~0.8 mm, averaged voxel
                Average_voxel_layers(4,:)  = mean(cortex{num_slice}{exp_case_num}(:,9:14),2);%L5,0.8~1.4 mm, averaged voxel
                Average_voxel_layers(5,:)  = mean(cortex{num_slice}{exp_case_num}(:,15:20),2);%L6,1.4~2.0 mm, averaged voxel
                if nCortical_Voxel > 20
                    Average_voxel_layers(6,:)  = mean(cortex{num_slice}{exp_case_num}(:,21:nCortical_Voxel),2);%WM,2.1~3.0 mm, averaged voxel
                end
            end
            for layer_idx = 1 : size(Average_voxel_layers,1)
                res_layer_line = abs(Average_voxel_layers(layer_idx,:));
                %         res_line2 = zscore(abs(squeeze(res_line(:,:,num_slice))),0,2);
                res_layer_line2 = res_layer_line - mean(res_layer_line,2);

                %         [s, f, t, power] = spectrogram(res_line2,hanning(128),96,1024,Fs,'yaxis','power');
                [s_voxel_line, f_voxel_line, t_voxel_line, line_layer_power] = spectrogram(res_layer_line2,hanning(128),96,1024,Fs,'yaxis','power');
                
                h=figure();
                add_offset = min(min(10*log10(abs(line_layer_power))))*(-1)/2;
                surf(t_line, f_line, 10*log10(abs(line_layer_power))+add_offset, 'EdgeColor', 'none');
                axis xy;
                axis tight;
                if Filter_On == 0
                    colorbar_max = 50;
                    colorbar_min = 0;
                else
                    colorbar_max = 100;
                    colorbar_min = -50;
                end
                colormap(jet);
                view(0,90);
                %     xlabel('Time(sec)');
                Xlim = get(gca, 'xlim');
                set(gca, 'XTick', linspace(Xlim(1), Xlim(end), 5));
                if TA == 640
                    xticklabels({'0s','160s','320s','480s','640s'})
                elseif TA == 768
                    xticklabels({'0s','192s','384s','576s','768s'})
                end
                colorbar;
                ylabel('Frequency(Hz)');
                set(gcf,'color','w');
                set(gca,'FontSize',font_size,'FontWeight','Bold');
                set(gca,'linewidth',axs_line_width);
                colormap(jet);
                c =colorbar('YTickLabel',{num2str(colorbar_min), num2str(colorbar_max)},'YTick', [colorbar_min colorbar_max]);
                c.Label.String = 'Power (dB)';
                caxis([colorbar_min colorbar_max])
                set(gcf, 'Position', [500, 500, 900, 500])
                %     strName3 = 'Blood Pressure Spectrogram (dB)';
                set(c,'YTickMode','manual');
                if layer_idx == 1
                    %             title('Layer #1');
                    strName4 = sprintf('Slice #%d, Line-scanning Layer #%d Spectrogram',num_slice,layer_idx);
                elseif layer_idx == 2
                    %             title('Layer #2-3');
                    strName4 = sprintf('Slice #%d, Line-scanning Layer #%d-%d Spectrogram',num_slice,layer_idx,layer_idx+1);
                elseif layer_idx > 2 && layer_idx < 6
                    %             title(sprintf('Layer #%d',layer_idx+1));
                    strName4 = sprintf('Slice #%d, Line-scanning Layer #%d Spectrogram',num_slice,layer_idx+1);
                elseif layer_idx == 6
                    %             title('WM');
                    strName4 = sprintf('Slice #%d, Line-scanning WM Spectrogram',num_slice);
                end
                
                title(strName4);
%                 strName40 = sprintf('Line-scanning #%d, Slice #%d, Layer #%d Spectrogram (dB)',str2num(data_num),num_slice,layer_idx);
                        strName40 = strcat('Line-scanning', '#', str2num(data_num),',',strName4);
                saveas(h,strName40,'png');
            end
        end
%% bugs exist!
%   h=figure();
%   nSingal = 2;
%   for index = 1 : nSingal
%       subplot (nSingal,1,index);
%       if index == 1
%           %for blood pressure
%           x1 = synchorized_blood_press(1,:);
%           [pxx_blood,f] = pwelch(x1,[],[],[],Fs1, 'onesided');
%           plot(f,10*log10(pxx_blood),'r-','LineWidth',plot_linewidth);
%           strName1 = sprintf('Blood Pressure PSD');
%           title(strName1);
%       elseif index == 2
%           %for respiratory press
%           x2 = synchorized_respiratory_freq(1,:);
%           [pxx_resp,f] = pwelch(x2,[],[],[],Fs1, 'onesided');
%           plot(f,10*log10(pxx_resp),'b-','LineWidth',plot_linewidth);
%           strName1 = sprintf('Respiratory Signal PSD');
%           title(strName1);
%       end
%       ylim([-60 60])
%       ylabel('PSD');
%       yticklabels('auto')
%       xticklabels('manual')
%       xlim([0 Fs1/2])
%       xticks(0:1:Fs1/2)
%       xticklabels({'0Hz','1Hz','2Hz','3Hz','4Hz','5Hz','6Hz','7Hz','8Hz','9Hz','10Hz'})
%       %     xlim([0 640*(1/TR)])
%       %     xticks(0:160*(1/TR):640*(1/TR))
%       %     xticklabels({'0','160 sec','320 sec','480 sec','640 sec'})
%       %     xlabel('Time(sec)');
%       set(gcf,'color','w');
%       set(gca,'FontSize',font_size,'FontWeight','Bold');
%       set(gca,'linewidth',axs_line_width);
%       set(gcf, 'Position', [500, 500, 900, 700])
%       
%       strName4 = sprintf('Blood and Resp #%d, PSD plots',str2num(data_num));
%       saveas(h,strName4,'png');
%   end
% end
%% xcorr vs. fMRI
if RS == 0
%         xc_lag = 100; %first trial
%         xc_lag = 200;
        xc_lag = 400;
else
        xc_lag = 500;
end

disp('x-corr cardiac vs. fmri')
    
for is = 1:Total_Slice
%     for ir = 1:length(Exp_Number_Array1)
        for iL = 1: size(cortex{is}{exp_case_num}, 2) %voxel size if 50um, 40 voxels, 100um, 20 voxels
            
            % for cardiac
            [xc_67{is}(exp_case_num, :, iL), lag_cardiac] ...
                = xcorr(zscore(cardiac_freq_67{exp_case_num}(:,1:TA*(1/TR)),0,2), zscore(cortex{is}{exp_case_num}(:, iL)), xc_lag, 'coeff');
%                 = xcorr(zscore(cardiac_freq_69{ir}), zscore(cortex{is}{ir}(:, iL)), xc_lag, 'coeff');

            % for respiratory
            [xc_01{is}(exp_case_num, :, iL), lag_respiratory] ...
                = xcorr(zscore(respiratory_freq_01{exp_case_num}(:,1:TA*(1/TR)),0,2), zscore(cortex{is}{exp_case_num}(:, iL)), xc_lag, 'coeff');

            
            %             [xc_05{is}(ir, :, iL), lag] = xcorr(zscore(ca_freq_05{ir}), zscore(ctx{is}{ir}(:, iL)), xc_lag, 'coeff');
            %             [xc_01{is}(ir, :, iL), lag] = xcorr(zscore(ca_freq_01{ir}), zscore(ctx{is}{ir}(:, iL)), xc_lag, 'coeff');
        end
%     end
    

    
    cc_67{is} = (xc_67{is}(exp_case_num, lag_cardiac==0, :));
    cc_01{is} = (xc_01{is}(exp_case_num, lag_respiratory==0, :));
    %     cc_05{is} = squeeze(xc_05{is}(:, lag==0, :));
    %     cc_01{is} = squeeze(xc_01{is}(:, lag==0, :));
    
    [~, im_67{is}] = max(xc_67{is}, [], 2);
    [~, im_01{is}] = max(xc_01{is}, [], 2);
    %     [~, im_05{is}] = max(xc_05{is}, [], 2);
    %     [~, im_01{is}] = max(xc_01{is}, [], 2);
end

%% corr across voxels and slices for cardiac
% for ir = 1:length(Exp_Number_Array1)
if Total_Slice  == 3
    h = figure; set(h, 'WindowStyle', 'docked', 'NumberTitle', 'off', 'Name', ['cardiac cc 67 r', num2str(Exp_Number_Array1(exp_case_num))])
    total_matrix = [squeeze(cc_67{1})'; squeeze(cc_67{2})'; squeeze(cc_67{3})'];
%     imagesc(1:Total_Slice, 1:length(cc_69{1}), total_matrix), colorbar, colormap(jet)
    imagesc(1:Total_Slice, 1:length(cc_67{1}), total_matrix'), colorbar, colormap(jet)
    %     imagesc(1, 1:length(im), squeeze(im)'), colorbar, colormap(jet)
% end


%% corr across voxels and slices for respiratory
% for ir = 1:length(Exp_Number_Array1)
    h = figure; set(h, 'WindowStyle', 'docked', 'NumberTitle', 'off', 'Name', ['resp cc 01 r', num2str(Exp_Number_Array1(exp_case_num))])
    total_matrix = [squeeze(cc_01{1})'; squeeze(cc_01{2})'; squeeze(cc_01{3})'];
%     imagesc(1:Total_Slice, 1:length(cc_69{1}), total_matrix), colorbar, colormap(jet)
    imagesc(1:Total_Slice, 1:length(cc_01{1}), total_matrix'), colorbar, colormap(jet)
    %     imagesc(1, 1:length(im), squeeze(im)'), colorbar, colormap(jet)
% end
end
 %%
 nCortical_Voxel = size(xc_67{is}, 3);
 Spatial_Resolution = 0.05;
% for ir = 1:length(Exp_Number_Array1)
    % for cardiac 
     h = figure; 
%      set(h, 'WindowStyle', 'docked', 'NumberTitle', 'off', 'Name', ['lag 2D 69 r', num2str(Exp_Number_Array1(ir))])
     for is=1:Total_Slice
         subplot(1,Total_Slice,is)
         imagesc(lag_cardiac*TR, 1:size(xc_67{is}, 3), squeeze(xc_67{is}(exp_case_num,:,:))')
         cmin = min(min(squeeze(xc_67{is}(exp_case_num,:,:))'));
         cmax = max(max(squeeze(xc_67{is}(exp_case_num,:,:))'));

         colormap(jet)
         colorbar
%          caxis([cmin cmax])
         caxis([-0.3 0.4])
%          caxis([-0.05 0.05])
         
         strTitle = sprintf('Slice #%d',is);
%          title(strTitle);
         set(gca,'FontSize',font_size,'FontWeight','Bold');
         set(gca,'linewidth',axs_line_width)
         array_index3 = 0:ceil(0.5/Spatial_Resolution):nCortical_Voxel;
         array_index3(1) = 1;
         yticks(array_index3)
         yticklabels({'0','0.5','1.0','1.5','2.0'});
         if is == 1
%             ylabel('Cortical depth (mm)')
         end
         if is == 2
         xlabel('Lag time (sec)')
         end
     end
     
 
%      for is=1:Total_Slice
%          subplot(1,3,is)
%          
%      end
     set(gcf,'color','w');
     if Total_Slice == 3
        set(gcf, 'Position', [500, 500, 900, 600])
     else
         set(gcf, 'Position', [500, 500, 400, 600])
     end
     strName0 = 'Cross correlation: Line-scanning fMRI vs. Cardiac Freq';
%      suptitle(strName0,'FontSize',font_size,'FontWeight','Bold');
     strName1 = sprintf('Cross correlation #%d, Line-scanning fMRI vs Cardiac Freq',Exp_Number_Array1(exp_case_num));
     saveas(h,strName1,'png');
     %%
     % for respiratory
      h = figure; 
%      set(h, 'WindowStyle', 'docked', 'NumberTitle', 'off', 'Name', ['lag 2D 69 r', num2str(Exp_Number_Array1(ir))])
     for is=1:Total_Slice
         subplot(1,Total_Slice,is)
         imagesc(lag_respiratory*TR, 1:size(xc_01{is}, 3), squeeze(xc_01{is}(exp_case_num,:,:))')
         cmin = min(min(squeeze(xc_01{is}(exp_case_num,:,:))'));
         cmax = max(max(squeeze(xc_01{is}(exp_case_num,:,:))'));

         colormap(jet)
         colorbar
         caxis([cmin cmax])
         
         strTitle = sprintf('Slice #%d',is);
         title(strTitle);
         set(gca,'FontSize',font_size,'FontWeight','Bold');
         set(gca,'linewidth',axs_line_width)
         array_index3 = 0:ceil(0.5/Spatial_Resolution):nCortical_Voxel;
         array_index3(1) = 1;
         yticks(array_index3)
         yticklabels({'0','0.5','1.0','1.5','2.0'});
         if is == 1
            ylabel('Cortical depth (mm)')
         end
         if is == 2
         xlabel('Lag time (sec)')
         end
     end
     
 
%      for is=1:Total_Slice
%          subplot(1,3,is)
%          
%      end
         set(gcf,'color','w');
     if Total_Slice == 3
        set(gcf, 'Position', [500, 500, 900, 600])
     else
         set(gcf, 'Position', [500, 500, 600, 600])
     end
     strName0 = 'Cross correlation: Line-scanning fMRI vs. Respiratory Freq';
     suptitle(strName0,'FontSize',font_size,'FontWeight','Bold');
     strName1 = sprintf('Cross correlation #%d, Line-scanning fMRI vs Respiratory Freq',Exp_Number_Array1(exp_case_num));
     saveas(h,strName1,'png');
end
%%


    
    
    












