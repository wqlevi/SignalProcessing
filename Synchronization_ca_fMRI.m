addpath(genpath('/Users/levi/Downloads/Single_Slice_LINE_data_for_paper/calcium_code'));
addpath('/Users/levi/Downloads/Single_Slice_LINE_data_for_paper');
% fMRI acquisition parameters
TR       = 1.5; %TR used in the EPI scans (seconds)
nTRs     = 410; %number of TRs
scanTime = nTRs*TR;
prestim  = 10; % seconds (calculate from TRs set in fmri parameters)

% figure properties
ftsz     = 10; %fontsize
lnw      = 2;
% import raw data
load('Scan17_TK.mat')
ca_raw = channels{1, 7}.data;
ca_raw = -ca_raw;
sti_raw = channels{1, 4}.data;
ca_sr = channels{1, 7}.samples_per_second;
sti_sr = channels{1, 4}.samples_per_second;
ca_t = round(0:1/ca_sr:(length(ca_raw)-1)/ca_sr);
sti_t = round(0:1/sti_sr:(length(sti_raw)-1)/sti_sr);

sti_index = find(sti_raw>1);
sti_trigger_idx = sti_index(1);
b4trigger_t = sti_trigger_idx/sti_sr; 
both_start_t = b4trigger_t - prestim/TR;

ca_newstart = ca_raw(both_start_t*ca_sr:(both_start_t+scanTime)*ca_sr); % now calcium matchs BOLD after crop

[P,Q] = rat(1/TR/ca_sr);
ca_newstart_sr = resample(ca_newstart,P,Q);

f_1 = figure();
re_t = (0:TR:(length(ca_newstart_sr)-1)*TR);
subplot(2,1,1)
plot(ca_t,ca_raw);
xlim([0 668])
xlabel('Time(sec)');
ylabel('amplitude(volts)');
title('Calcium raw data');
grid on;

subplot(2,1,2)
plot(re_t,ca_newstart_sr);% resampled and croped calcium (seconds)
xlim([1 length(ca_newstart_sr)-1])
xlabel('Time(sec)');
ylabel('amplitude(volts)');
title('Calcium matched with fMRI');
grid on;
saveas(f_1,'raw plot','png');

var        =  ca_newstart;
varname    = 'adjusted Ca';
tDuration  = nTRs*TR; %seconds
fs         = ca_sr;
windowSize = 2;  % s  
[tf freqs times] = timefreq(var,fs,'wletmethod','dftfilt3','winsize',windowSize*fs,'ntimesout',round(tDuration/windowSize*2),'freqs',[0,20]);
 figure;      % subplot (311);
        tf = abs(tf); % get the power
        times = times/1000; %mse --> sec
        imagesc(times,freqs,tf);
        set(gca,'ydir','normal','xlim',[times(1) times(end)],'ylim',[0.1 20]);
%         xlabel(sprintf('Time in sec')); 
        ylabel('Frequency in Hz'); grid on;
        set(gca,'XTick',[0:60:60*tDuration]);
%         set(gca,'XTickLabel',[0:60:60*tDuration]);
        title(strcat('calcium spectrum'));
        %title(strcat(num2str(minNum(1,iMin)),'min after surgery'));  %%%%%
        %for automation use the iMin
        h=colorbar;
        ylabel(h,'Power');
        caxis([0, 30]); 
