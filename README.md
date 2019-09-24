# SignalProcessing
CALCIUM
This repo contains signal processing with calcium/BOLD signals

Evoked:
- Reverse raw data in time course(-1*raw data);
- Eliminate freq <1 Hz;
- Focus on 3-6Hz(Responses, 6Hz has cardiac overlapped w/ harmonics);

Resting-State:

fMRI
Normally, BOLD signal is observed via either AFNI(from NIH) or try to analysis on Matlab. Now one more choice is to do so using Nipy(http://pysurfer.github.io/auto_examples/index.html)

From Nipy, signals can be ploted on a 3D reconstruction model of human brain.(Which indicates animal ones)
