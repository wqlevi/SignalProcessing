% multi-animal date averaging
% date = ['2809';'0510';'2707'];%opto
% date = ['1407';'3107'];
% date = ['1407';'0310';'0910'];% sole side
% date = ['1407';'0310';'0910';'1607'];% sole side tk & Hang's
% date = ['1407';'0310';'0910';'1607'];% rs & Hang's
date = ['2707';'2809';'0510';'1008'];% opto & Hang's
% date = ['2707';'1008'];
% date = ['0910';'2707'];
% date = ['3107';'2707'];
nanimals = size(date,1);
nslice = 2;
TR = 0.1;
nY = 32;
ndur = 200;
t_fmri = 0:TR:(nY*ndur-1)*TR;
t_epoch = [0:TR:(ndur-1)*TR];
Font = 15;
layer_idx = [1,2,2,2,2,3,3,4,4,4,4,4,4,5,5,5,5,5,5,5];
layer_name = {'L1','L2/3','L4','L5','L6'};
nvoxels = 20;
% load data here
for ia = 1:nanimals
    data_path = ['D:\',date(ia,:),'2020_yi\figs'];
    cd(data_path);
    tmp = load('trial_ave.mat');
    fmri_tc{ia} = tmp.fmri_time_temp;
    tmp1 = load('epoch_ave.mat');
    epoch_tc{ia} = tmp1.fmri_cm_epoch;
    tmp2 = load('norm_fmri.mat');
    raw_fmri{ia} = tmp2.norm_cm;
%     tmp3 = load('pxx.mat');
%     pxx{ia} = tmp3.pxx;
    tmp4 = load('f.mat');
    f = tmp4.f;
%     tmp5 = load('clm.mat');
%     clm{ia} = tmp5.norm_cortical_map;
    tmp6 = load('snr_sp.mat');
    snr_sp{ia} = tmp6.tsnr_temp1;
    tmp7 = load('all_snr.mat');
    all_snr{ia} = tmp7.tsnr_temp;
    
    tmp8 = load('cort_idx.mat');
    idx{ia} = tmp8.half_idx;

end
for ia = 1:nanimals
    for is = 1:nslice
        for ir = 1:size(raw_fmri{ia},2)
            all_raw(ia,is,:,:) = raw_fmri{ia}{1}{is};
        end
        ave_tmp = mean(fmri_tc{ia}{1}{is},1);
        all_tmp(ia,is,:,:) = fmri_tc{ia}{1}{is};
        std_tmp = std(fmri_tc{ia}{1}{is},1);
        ave(ia,is,:) = ave_tmp;
        ave1(ia,is,:,:) = epoch_tc{ia}{1}{is};
        stdv(ia,is,:) = std_tmp;
%         for ir = 1:size(pxx{ia},1)
%         ma = max(pxx{ia}(ir,is,0.01<f&f<0.5));
%         mi = min(pxx{ia}(ir,is,0.01<f&f<0.5));
%         pxx_norm{ia}(ir,is,:) = (pxx{ia}(ir,is,0.01<f&f<0.5)-mi)/(ma-mi);
%         end
%         fmri_map(ia,is,:,:) = clm{ia}(:,:,1,is);
    end
end
% assign data from cell to matrix
all_ave = squeeze(mean(ave,1));%all animal Percentage ave
ave_raw = squeeze(mean(all_raw,1));% all animal RAW ave
all_ave_ep = squeeze(mean(ave1,1));
all_std = squeeze(mean(stdv,1));
ave_tmp = squeeze(mean(all_tmp,1));
for ia = 1:nanimals
    pxx_norm_tmp(ia,1,:) = mean(pxx_norm{ia}(:,1,:),1);
end
pxx_norm_all = squeeze(mean(pxx_norm_tmp,1));
% ave_map = mean(fmri_map,1);

%% Time course
% timecourse of ave 
for is = 1:nslice  
    h = figure;
    %errorbar(t_fmri,100*all_ave(is,:),all_std(is,:),'Color',[.7,.7,.7],'LineWidth',1);
%     hold on;
    plot(t_fmri/60,100*all_ave(is,:),'k','LineWidth',1);
    xlim([t_fmri(1),t_fmri(end)/60]);
    xlabel('Time (min)');
    ylabel('BOLD(%)');
    set(gcf,'Position',[500 500 980 300]);
    ax = gca;
    ax.XAxis.FontSize = Font;
    ax.YAxis.FontSize = Font;
    ax.FontWeight = 'bold';
    ax.LineWidth = 2;
    ylim([-1,2]);
    box off;
end
Clr = flip(colormap(jet(20)),1);
for is = 1:nslice
    h = figure;
    for iv = size(ave1,3):-1:1
    plot(t_epoch,100*squeeze(all_ave_ep(is,iv,:)),'Color',Clr(iv,:));
    hold on;
    end
    ylabel('Percentage(%)');
    colormap(Clr);
    ax = gca;
    ax.XAxis.FontSize = Font;
    ax.YAxis.FontSize = Font;
    ax.FontWeight = 'bold';
    ax.LineWidth = 1;
    colorbar('Ticks',[0,1],'Direction','reverse','TickLabels',{'Surface','Bottom'});
    box off;
end
%% tSNR whole LS
% formatting snr 
for is = 1:nslice
    figure;
    for ia = nanimals:-1:1
        if ia ~= 1
            for ir = 1:size(all_snr{ia},2)
                x_id(ir,ia,:) = 1:128;
                x_id(ir,ia,:) = x_id(ir,ia,:)+(idx{1}(:,ir,is)-idx{ia}(:,ir,is));% align cortex ref. 1st animal
            end
            idx_min = min(min(x_id(:,2:end,end))); % R-edge of x-axis
            idx_max = max(max(x_id(:,2:end,1)));% L-edge of x-axis
            for ir = 1:size(all_snr{ia},2)
                all_snr_1{ia} = all_snr{ia}((x_id(ir,ia,:)>=idx_max & x_id(ir,ia,:)<=idx_min),:,:);% trim arrays for rest animals
                plot(squeeze(all_snr_1{ia}(:,ir,is)),'b-');
                hold on;
            end
        else 
            all_snr_1{ia} = all_snr{ia}(idx_max:idx_min,:,:);
            plot(squeeze(all_snr_1{ia}(:,ir,is)),'r-');
            hold on;
        end
    end
end

for is = 1:nslice
    for ia = 1:nanimals
        if ia == nanimals
            pre_snr_0(1,:,is) = squeeze(mean(all_snr_1{ia}(:,:,is),2));
        else 
            cur_snr_0(ia,:,is) = squeeze(mean(all_snr_1{ia}(:,:,is),2));
        end
    end
end

ave_cur_snr_0 = mean(cur_snr_0,1);
std_cur_snr_0 = std(cur_snr_0,0,1);
std_pre_snr_0 = std(all_snr_1{nanimals},0,2);

clr_map_line = [0,0,255;255,51,51]./255;
clr_map_shadow = [153,204,255;255,153,153]./255;


for is = 1:nslice
    figure; 
    for ia = 1:nanimals
        for ir = 1:size(all_snr{ia},2)
            hold on;
            errorbar(1:size(all_snr_1{ia}),squeeze(pre_snr_0(:,:,is)),squeeze(std_pre_snr_0(:,:,is)),'Color',clr_map_shadow(2,:),'LineWidth',1);
            hold on;
            errorbar(1:size(all_snr_1{ia}),squeeze(ave_cur_snr_0(:,:,is)),squeeze(std_cur_snr_0(:,:,is)),'Color',clr_map_shadow(1,:),'LineWidth',1);
%             xlim([1,size(all_snr_1{ia})]);
            xlabel('Cortical depth(voxel)');    
            ylabel('tSNR');
        end
        plot(1:size(all_snr_1{ia}),squeeze(ave_cur_snr_0(1,:,is)),'Color',clr_map_line(1,:),'LineWidth',3);
        plot(1:size(all_snr_1{ia}),squeeze(pre_snr_0(1,:,is)),'Color',clr_map_line(2,:),'LineWidth',3);
    end
    ax = gca;
    ax.XAxis.FontSize = Font;
    ax.YAxis.FontSize = Font;
    ax.FontWeight = 'bold';
    ax.LineWidth = 1;
end
%% tSNR

%ave animals
for is = 1:nslice
    for ia = 1:nanimals
        if ia == nanimals
            pre_snr(1,:,is) = squeeze(mean(snr_sp{ia}(:,:,is),2));
        else 
            cur_snr(ia,:,is) = squeeze(mean(snr_sp{ia}(:,:,is),2));
        end
    end
end

ave_cur_snr = mean(cur_snr,1); % average of current animals
std_cur_snr = std(cur_snr,1);
std_pre_snr = std(snr_sp{nanimals},0,2);


% separate trials
for is = 1:nslice
    figure; 
    for ia = 1:nanimals
%         plot(1:20,squeeze(snr{ia}(:,1,is)),'Color',clr_map_line(ia,:),'LineWidth',3);% ave of raw data trail
        for ir = 1:size(snr_sp{ia},2)
            hold on;
%             if ia == nanimals % for last animal(previous one)
%                 cl_i = 2;
%             
%             else
%                 cl_i = 1;
%             end
%             plot(1:20,squeeze(snr_sp{ia}(:,ir,is)),'Color',clr_map_shadow(cl_i,:));
            errorbar(1:20,squeeze(pre_snr(:,:,is)),squeeze(std_pre_snr(:,:,is)),'Color',clr_map_shadow(2,:),'LineWidth',1);
            hold on;
            errorbar(1:20,squeeze(ave_cur_snr(:,:,is)),squeeze(std_cur_snr(:,:,is)),'Color',clr_map_shadow(1,:),'LineWidth',1);
            xlim([1,20]);
            xlabel('Cortical depth(voxel)');    
            ylabel('tSNR');
        end
        plot(1:20,squeeze(ave_cur_snr(1,:,is)),'Color',clr_map_line(1,:),'LineWidth',3);
        plot(1:20,squeeze(pre_snr(1,:,is)),'Color',clr_map_line(2,:),'LineWidth',3);
    end
    ax = gca;
    ax.XAxis.FontSize = Font;
    ax.YAxis.FontSize = Font;
    ax.FontWeight = 'bold';
    ax.LineWidth = 1;
%     title(['tSNR comparison slice#',num2str(is)]);
end
%% 2D colormaps
for is = 1:nslice
    h = figure;
    imagesc(t_epoch,[],100*squeeze(all_ave_ep(is,:,:)));
    colormap jet;
        switch is
        case 1
        caxis([-0.2,3.5]);
        case 2
        caxis([-.2,0.4]);
    end
    h1 = colorbar;
    ylabel(h1,'Percentage (%)','FontSize',Font,'FontWeight','bold');
    ax = gca;
    ax.XAxis.FontSize = Font;
    ax.YAxis.FontSize = Font;
    ax.FontWeight = 'bold';
    ylabel('Cortical depth(voxels)');
    xlabel('time(s)');
end


% filtering & 10 percent pinpoint[EPOCH]
[b_1,a_1] = butter(1,0.12,'low');
len_ep = size(all_ave_ep,3);

for is = 1:nslice
    for iv  = 1:nvoxels
        fmri_ep_tmp{is}(iv,:) = all_ave_ep(is,iv,:);
        temp = [fliplr(fmri_ep_tmp{is}(iv,1:len_ep/4)),fmri_ep_tmp{is}(iv,:),fliplr(fmri_ep_tmp{is}(iv,end-len_ep/4+1:end))];
        temp_filt3 = filtfilt(b_1,a_1,temp);
        fmri_ep{is}(iv,:) = temp_filt3(:,len_ep/4+1:end-len_ep/4);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%test filtering
% test = fmri_tc{1}{1}{1};
Fs= 1/TR;

[b_2,a_2] = butter(1,[0.01,0.1]/((1/TR)/2),'bandpass');
order    = 4096; 
Fc1  = 0.01;   %low cut frequency in Hz 0.01 
Fc2 = 0.1;   %intend 0.1 Hz but consider transient band
N = order;
beta = 0.005;
win = kaiser(N+1, beta); %using kaiser window
flag = 'scale';  % Sampling Flag

% Calculate the coefficients using the FIR1 function.
b  = fir1(N, [Fc1, Fc2]/(Fs/2), 'bandpass', win, flag); % causing dampened ampitude and phase shift!!

fvtool(b,1,'Fs',Fs)
Hd = dfilt.dffir(b);
% MAINLY USED
for is = 1:nslice
    test = squeeze(ave_raw(is,:,:));
    test3 = test ./ max(max(test));
    % test3 = squeeze(test);
    for iv  = 1:nvoxels
    temp2 = [fliplr(test3(iv,2:order+1)),test3(iv,:),fliplr(test3(iv,end-order:end-1))];
    temp1 = [fliplr(test3(iv,1:(6400/4))),test3(iv,:),fliplr(test3(iv,end-(6400/4)+1:end))];
    temp_filt2 = filtfilt(b_2,a_2,temp1);
    test1(iv,:) = temp_filt2(:,6400/4+1:end-6400/4);
    temp_filt3 = filter(Hd,temp2);
    delay = mean(grpdelay(Hd));
    temp_filt4 = circshift(temp_filt3,(-1)*delay,2);
    test2(iv,:) = temp_filt4(:,order+1:end-order);
    end
    % plot(t_fmri,squeeze(mean(fmri_tc{1}{1}{1},1)),'k',t_fmri,squeeze(mean(test1(1,:),1)),'r',t_fmri,squeeze(mean(test2(1,:),1)),'b');
    % legend({'Percentage','IIR raw','FIR raw'});
    h = figure;
    imagesc(t_fmri/60,[],test1 ./ max(max(test1)));colormap jet;colorbar;
    ylabel('IIR');
    caxis([0,0.8]);
    set(gcf,'Position',[500,500,900,300]);
    ax = gca;
    ax.XAxis.FontSize = Font;
    ax.YAxis.FontSize = Font;
    ax.FontWeight = 'bold';
    ylabel('Cortical depth(voxels)');
end
% AVE FROM PREVIOUS DATA
% for is = 1:nslice
% h = figure;
% imagesc(t_fmri,[],squeeze(ave_map(1,is,:,:)));colormap jet;colorbar;
% ylabel('IIR');
% caxis([0.5,0.8]);
% set(gcf,'Position',[500,500,900,300]);
% ax = gca;
% ax.XAxis.FontSize = Font;
% ax.YAxis.FontSize = Font;
% ax.FontWeight = 'bold';
% ylabel('Cortical depth(voxels)');
% end
subplot(224)
imagesc(t_fmri,[],test2 ./ max(max(test2)));colormap jet;colorbar;
ylabel('FIR');
caxis([0.45,01])
set(gcf,'Position',[500,500,900,300]);
%%%%%%%%%%%%%%%%%
for is = 1:nslice
    h = figure;
    for il = 1:5
    subplot(5,1,il);
    h1=plot(t_epoch,100*squeeze(mean(all_ave_ep(is,layer_idx==il,:),2)),'k');
    hold on;
    [~,pks1] = find(mean(fmri_ep{is}(layer_idx == il,:),1)>0.5*max(mean(fmri_ep{is}(layer_idx == il,:),1)),1);
%     [~,pks2] = find(mean(fmri_ep{is}(layer_idx == il,:),1)>0.1*max(mean(fmri_ep{is}(layer_idx == il,:),1)),1);% 10% max
%     plot(pks1/10,100*mean(fmri_ep{is}(layer_idx == il,pks1),1),'r*','MarkerSize',10);% 50% max
%     plot(pks2/10,100*mean(fmri_ep{is}(layer_idx == il,pks2),1),'g*','MarkerSize',5);% 10% max
    pks_idx(is,il) = pks1;
    ax = gca;
    if il ~=5 
    ax.XAxis.Visible = 'off';
    end
    box off;
    if il == 1
    set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%     legend({'50% peak amplitude'});
%     legend('boxoff')
    end
    ax.XAxis.FontSize = Font;
    ax.YAxis.FontSize = Font;
    ax.LineWidth = 1;
    ax.FontWeight = 'bold';
    set(gcf,'Position',[490 86 450 670]);
    xlabel('time(s)');
    ylabel(layer_name{il});
    ytickformat('%g%%');
    max_ep(is,il) = max(mean(fmri_ep{is}(layer_idx == il,:),1));
    min_ep(is,il) = min(mean(fmri_ep{is}(layer_idx == il,:),1));
    end
end

for is = 1:nslice
    h = figure;
    bar([100*max_ep(is,:);100*min_ep(is,:)]');
    ylabel('Percentage(%)');
%     ylim([-2,2])
    xticklabels({'L1','L2/3','L4','L5','L6'});
    title(['Amplitude of Positive & Negative BOLD response slice#',num2str(is)]);
    saveas(h,['Amplitude of Positive & Negative BOLD response slice#',num2str(is),'.jpg']);
end

% separated bars
for is = 1:nslice
    h = figure;
    bar(1:5,pks_idx(is,:)/10);
    ax = gca;
    ax.XTickLabel = layer_name;
    ax.YLabel.String = 'On-set time(sec,50% amplitude)';
    box off;
end
% all bars
h = figure;
bar(pks_idx'/10);
ylim([0,5]);
ax = gca;
ax.XTickLabel = layer_name;
ax.YLabel.String = 'On-set time(sec,50% amplitude)';
box off;
legend({'Slice1','Slice2'},'Location','northeastoutside');
legend boxoff;
