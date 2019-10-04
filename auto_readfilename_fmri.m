
% Path to data and utilities
path_ca = 'C:\Users\wangqi\Documents\Lab\Data\05302019_Hang\05302019_ca';
addpath ('C:\Users\wangqi\Documents\Lab\Demo\Calcium');
% AutoRead files
%  for i = 1:5
%      i_str = num2str(i);
%      name(:,i) = strcat('Scan_',i_str,'_TK.mat');
%  end
cd (path_ca)
list_ca = ls("*Scan*.acq");
for j = 1:size(list_ca,1) %divide filename
    number_name(:,j) = split(list_ca(j,:),'_');
    string_ca(j) = convertCharsToStrings(number_name{3,j}(:));
    strcmp(mode_paradigm, 'rs')
%     result_str(j) = contains(string_ca(j),'TK'); % if it's TK scan
%     if result_str(j)==1
%         mode(:,j) = "Task";
%     else
%         mode(:,j) = "Resting-State";
%     end
 end
 
%
% Manually pick data in user interface
%     [f,p] = uigetfile([{'fid*;ser*', ...
%             'Bruker raw data files';'*.*','All Files'}], ... 
%             'Select fid-file for reading')
%     if f == 0
%         arr = -1;
%         return;
%     end
%     path = strcat(p,f);
% 
% %
% Experimrntal setting
if TR50 == 1
    TR = 0.05; %s
    b4trig = 20; %prestim numbers
    nimg = 400;
%     total_scan_time = 640; %second
else
    TR = 0.1; %s
    b4trig = 10; %prestim numbers
    nimg = 200;% numbers of img per epoch
%     total_scan_time = 640; %second
end
%

nY = 32;% number of epoches
total_scan_time = nimg*TR*nY;% 640sec per scan

fmri_tmp = zero()
