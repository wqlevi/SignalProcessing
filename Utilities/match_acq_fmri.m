function [acqData, fmri, begin, fin] = match_acq_fmri(acqData, fmri, TR, b4trig)
% match the calcium and fmri time courses. fMRI starts 'b4trig' TR before the onset of the first trigger (stim channel)
% get b4trig using the shell script 'list_preBaselineNum_nTR_TR'

nOmit = 0;
%     fmri = fmri(1+nOmit:end);

    for i=1:length(acqData.channels) % find Stimulator Output  
        if strcmp('Stimulator Output',acqData.channels{i}.name)
            i_trig = find( abs(acqData.channels{i}.data) > 0.5, 1); % find start point with value >0.5
            t_trig = 0:1/acqData.channels{i}.samples_per_second:(length(acqData.channels{i}.data)-1)/acqData.channels{i}.samples_per_second;
        end
    end

    fmri_longer = 0;
    for i=1:length(acqData.channels)
        fs = acqData.channels{i}.samples_per_second;
        t_chan = 0:1/fs:(length(acqData.channels{i}.data)-1)/fs;
        i_begin = find(t_chan>=t_trig(i_trig),1);
        i_begin = round(i_begin - fs*b4trig*TR + fs*nOmit*TR); % if fs isn't an integer - round (Jans data)
        begin(i) = i_begin;
        if max(t_chan)-t_chan(i_begin)+1>(length(fmri))*TR
            i_end = round(i_begin - 1 + fs*(length(fmri)-nOmit)*TR);
            fin(i) = i_end;
        else
            i_end = length(t_chan);
            fin(i) = i_end;
            fmri_longer = 1;
        end
        acqData.channels{i}.data = acqData.channels{i}.data(i_begin:i_end);
    end

    if fmri_longer
        t_fmri = 0:TR:(length(fmri)-1)*TR;
        fmri = fmri(t_fmri<=max(t_chan)-t_chan(i_begin));
    end
end
