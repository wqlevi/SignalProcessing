data_path = 'D:\RS_ca\';
data_str = {'01222019','03192019','03242019','05232019','08092019','09122019','10062018','10092018','10252018','11242019'};
suf_str = 'ca_RS\';
for it = 1:size(data_str,2)
switch it
    case 1
        scans{it} = [26,27,28,29,30,31,32,33,35,36,38];% 01222019
    case 2
        scans{it} = [15,16];% 03192019
    case 3
        scans{it} = [35,37,39,40,41,44];% 03242019
    case 4
        scans{it} = [17,19,21];% 05232019
    case 5
        scans{it} = [13,15,19,22,24,31];% 08092019
    case 6
        scans{it} = [13,15,17,19,21,23,25,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,68,70,72,74,76,78];% 02022019
    case 7
        scans{it} = [30,36,38,40,42,43,45,47,51,53,55,57,59,61,63,73];% 10062018
    case 8
        scans{it} = [22,26];% 10092018
    case 9
        scans{it} = [41,42,43,44,47,48,51,52,53,54,55,56];% 10252018
    case 10
        scans{it} = [11,21,23,25,28,30];% 11242019
    end
        
end     
        












TR = 0.1;
fmri_duration = 640;
fmri_dummy = ones(fmri_duration/TR,1);
prestim = 10;

for it = 1:size(data_str,2)
for ir  = 1:length(scans{it})
    ca = load([data_path,data_str{it},suf_str,'scan_',num2str(scans{it}(:,ir)),'.mat']);
    
    [data_match,fmri_dummy,beg,fin] = match_acq_fmri(ca,fmri_dummy,TR,prestim);
    fs_ca{it} = data_match.channels{7}.samples_per_second;
    ca_match{it}{ir} = -data_match.channels{7}.data';
    
end
end
t_ca = [0:1/fs_ca{it}:(length(ca_match{it}{ir})-1)/fs_ca{it}];
for it = 1:size(data_str,2)
    for ir = 1:length(scans{it})
          ca_demean{it}{ir} = ca_match{it}{ir} - mean(ca_match{it}{ir});
    end
end
win = fs_ca{1}*2;
for it = 1:size(data_str,2)
    for ir = 1:length(scans{it})
    % for ir = trails
        [pxx(:,ir),f] = pwelch(ca_demean{it}{ir},length(ca_demean{it}{ir})/4,[],[],fs_ca{1}); 
        pxx_temp(:,ir) = zscore(10*log10(pxx(f>0.01&f<10,ir)),[],1);
        pxx_n{it}(:,ir) = pxx_temp(:,ir) - min(pxx_temp(:,ir),[],1);
    end
end
for it = 1:size(data_str,2)
    f_n = f(f>0.01&f<10);
    pxx_ave  = mean(pxx_n{it},2);   
    pxx_min = min(pxx_n{it},[],2);
    pxx_max = max(pxx_n{it},[],2);

    cd([data_path,data_str{it},suf_str])
    h = figure;
    k = fill([f_n;flipud(f_n)],[pxx_min;flipud(pxx_max)],[.9 .9 .9],'linestyle','none');
    line(f_n,pxx_ave,'Color','black','LineWidth',2);
    legend('all trails','mean of all trails');
    ylim([0,15]);
    xlabel('Frequency(Hz)');
    ylabel('Power(a.u.)');
    title('Averaged Calcium PSD');
    box off;
    text(7,12,['animal# ',data_str{it}],'FontWeight','bold');
    saveas(h,['PSD for animal#',data_str{it},'.png']);
    clear pxx pxx_temp 
end
 
%ave_all
pxxx_all = cell2mat(pxx_n);
clear pxx_n
pxxx_ave = mean(pxxx_all,2);
pxxx_min = min(pxxx_all,[],2);
pxxx_max = max(pxxx_all,[],2);
h = figure;
k = fill([f_n;flipud(f_n)],[pxxx_min;flipud(pxxx_max)],[.9 .9 .9],'linestyle','none');
line(f_n,pxxx_ave,'Color','black','LineWidth',2);
legend('all trails','mean of all trails');
ylim([0,12]);
xlabel('Frequency(Hz)');
ylabel('Power(a.u.)');
title('Averaged Calcium PSD');
box off;
text(7,12,['animal# ',data_str{it}],'FontWeight','bold');
saveas(h,['PSD for all animals']);