% multi-animal date averaging
date = ['2809';'0510';'2707'];%opto
% date = ['1407';'3107'];
% date = ['1407';'0310';'0910'];% sole side
% date = ['1407';'0310';'0910';'1607'];% sole side tk & Hang's
% date = ['1407';'0310';'0910';'1607'];% rs & Hang's
% date = ['2707';'2809';'0510';'1008'];% opto & Hang's
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
    end
end
% assign data from cell to matrix
all_ave = squeeze(mean(ave,1));%all animal Percentage ave
ave_raw = squeeze(mean(all_raw,1));% all animal RAW ave
all_ave_ep = squeeze(mean(ave1,1));
all_std = squeeze(mean(stdv,1));
ave_tmp = squeeze(mean(all_tmp,1));
%% RS tsnr whole
clr_map_line = [0,0,255;255,51,51]./255;
clr_map_shadow = [153,204,255;255,153,153]./255;

for is = 1:nslice
    for ia = nanimals:-1:1
            for ir = 1:size(all_snr{ia},2) 
                x_id{ia}(ir,is,:) = 1:128;
                x_id{ia}(ir,is,:) = x_id{ia}(ir,is,:)+(idx{1}(:,1,is)-idx{ia}(:,ir,is));% re-aligned cortex ref. 1st animal
            end
    end
end
idx_min = 119;
idx_max = 18;
for is = 1
    figure;
    for ia = nanimals:-1:1
            for ir = 1:size(all_snr{ia},2)
                all_snr_1{ia}(:,ir,is) = all_snr{ia}((x_id{ia}(ir,is,:)>=idx_max & x_id{ia}(ir,is,:)<=idx_min),ir,is);% trim arrays for rest animals
                if ia == 4
                    plot(squeeze(all_snr_1{ia}(:,ir,is)),'r-');
                    hold on;
                else 
                    plot(squeeze(all_snr_1{ia}(:,ir,is)),'b-');
                end
            end 
    end
end


for is = 1
    for ia = 1:nanimals
        if ia == nanimals
            pre_snr_1(1,:,is) = squeeze(mean(all_snr_1{ia}(:,:,is),2)); % ave of pre ani
        else 
            cur_snr_1(ia,:,is) = squeeze(mean(all_snr_1{ia}(:,:,is),2));
            disp(num2str(ia))
        end
    end
end

ave_cur_snr_1 = mean(cur_snr_1,1); % average of current animals
std_cur_snr_1 = std(cur_snr_1,1);
std_pre_snr_1 = std(all_snr_1{nanimals},0,2);

for is = 1
    figure;
    errorbar(1:size(pre_snr_1(1,:,is),2),squeeze(pre_snr_1(1,:,is)),squeeze(std_pre_snr_1(:,:,is)),'Color',clr_map_shadow(2,:),'LineWidth',1);
    hold on;
    xlim([1,size(pre_snr_1(1,:,is),2)]);
    ylim([0,70]);
    errorbar(1:size(pre_snr_1(1,:,is),2),squeeze(ave_cur_snr_1(:,:,is)),squeeze(std_cur_snr_1(:,:,is)),'Color',clr_map_shadow(1,:),'LineWidth',1);
    plot(1:size(pre_snr_1(1,:,is),2),squeeze(pre_snr_1(1,:,is)),'Color',clr_map_line(2,:),'LineWidth',3);
    hold on;
    plot(1:size(pre_snr_1(1,:,is),2),squeeze(ave_cur_snr_1(:,:,is)),'Color',clr_map_line(1,:),'LineWidth',3);
    set(gca,'LineWidth',2);
    box off;
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


for is = 1:nslice
    figure;
    for ia = 1:nanimals
        for ir = 1:size(all_snr{ia},2)
            if ia == nanimals 
                plot(squeeze(snr_sp{ia}(:,ir,is)),'r-');
            else
                plot(squeeze(snr_sp{ia}(:,ir,is)),'b-');
            end
            hold on;
        end
    end
end

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
    ax.LineWidth = 2;
%     title(['tSNR comparison slice#',num2str(is)]);
end