function [outsig] = ReplaceNan2Zero(insig)
% The function ReplaceNan2Zero will transform NaN elements in array into
% numerical zero.

tn = isnan(insig);
tz = isempty(tn);
% tf = isinf(insig);
% tff = find(tf);
    if ~tz
        disp('There exists NaN in array');
        insig(isnan(insig)) = 0;
        outsig = insig;
    else
        dis('There exists no NaN in array');
    end
end