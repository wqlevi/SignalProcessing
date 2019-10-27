function fid = BrReadFid(path)

% Reads Bruker raw data file.
% Reads the fid-file path without reading or expecting any additional
% information. The size of the resulting row vector depends only on the
% size of the file.

% Written by Rolf Pohmann
% Max Planck Institute for Biological Cybernetics, Tübingen, Germany
%###################Qi's editing#############################
% Open a window and show path of file, according to its name in format
%path
% Nargin is number of input parameters.

if nargin < 1
    [f,p] = uigetfile([{'fid*;ser*;2dsqe*', ...
            'Bruker raw data files';'*.*','All Files'}], ... 
            'Select fid-file for reading');
    if f == 0
        arr = -1;
        disp('No data selected');
        return;
    end
% Connecting two strings p and f
    path = strcat(p,f);
end
%% Open file and start reading
file = fopen(path,'r');
if file == -1
    arr = -2;
    return
end
% fread(,'int32') read binary data into int32 precision
fid = fread(file,'int32');
fclose(file);
% reshape allow reconstruct data in fid with format of 2 row and whatever
% columns
fid = reshape(fid,2,[]);
fid = complex(fid(1,:),fid(2,:));
fiid = fft(fid);
plot(fiid);
% fiid = abs(fiid);
% loglog(x,fiid);
%return;
end