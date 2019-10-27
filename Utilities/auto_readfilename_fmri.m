function [TK_test,RS_test] = auto_readfilename_fmri(path_ca)
% Path to data and utilities
% path_ca = 'C:\Users\wangqi\Documents\Lab\Data\05302019_Hang\05302019_ca';

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%   Usage: Input 'path_ca'                          %
%          Output trails number of 'RS' and 'TK'    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% AutoRead files
cd (path_ca)
list_ca = ls("*Scan*.acq");
TK_numbers = [];
RS_numbers = [];
TK_test = [];
RS_test = [];
for j = 1:size(list_ca,1) %divide filename
    file_name(:,j) = split(list_ca(j,:),'_');
    mode_name(j,:) = file_name{3,j};
    index_total(j,:) = convertCharsToStrings(file_name{2,j});
    if strfind(mode_name(j,:),'TK')
        TK_numbers = [TK_numbers,index_total(j,:)];
    elseif strfind(mode_name(j,:),'RS')
        RS_numbers = [RS_numbers,index_total(j,:)];
    end   

 end
    for i = 1:size(TK_numbers,2)
    TK_test = [TK_test,str2num(TK_numbers(:,i))];
    end
    for k = 1:size(RS_numbers,2)
    RS_test = [RS_test,str2num(RS_numbers(:,k))];
    end
    
 end