# SignalProcessing
##CALCIUM
This repo contains signal processing with calcium/BOLD signals

Evoked:
>- Reverse raw data in time course(-1*raw data);
>- Demean data, to eliminate 0Hz strong signal, and maybe a 10*log10(x) to be used to minimized dramatic signal drop; 
>- FFT timecourse signal into spectrum;
>- Match to fmri by cut off signal before stimulation
>- Extract targeted frequency band from spectrum, used as X-correlation input;
>- Stimulation freq = 3.33Hz, and 1Hz for respiratory, 6Hz for cardiac; (harmonics exists)

Resting-State:
>- Reverse raw data in time course;
>- *Demean* data, and plot an *pwelch* spectrum to visualize the frequency-power profile.(*10*log10(x) applied);
>- FFT timecourse into time-frequency decomposition(*Power spectrogram*);
>- Extract frequency band with salient power;
>- Apply an low BPF to power profile and see the correlation with BOLD;

##fMRI
Normally, BOLD signal is observed via either AFNI(from NIH) or visualized on Matlab(Brucker packages). Now one more choice is to do so using Nipy(http://pysurfer.github.io/auto_examples/index.html)

From Nipy, signals can be ploted on a 3D reconstruction model of human brain.(Which indicates animal ones)
