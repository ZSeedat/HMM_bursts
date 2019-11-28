# HMM_bursts
MATLAB code used for burst detection in source-localised MEG data.
This code is dependent on the HMM Toolbox created by OHBA. Go to :https://github.com/OHBA-analysis/HMM-MAR to download the toolbox.
There is also a detailed wiki explaining how to use each function.

# Scripts
- hmm_burst_detect.m finds bursts in source-localised parcellated MEG data.
- method_comparison_plot.m selects the burst state from the output of the HMM and plots it along with its spectrum. It also plots the bursts found using a simple thresholding technique for comparison.
- hmm_burst_features.m calculates the average burst count, duration, amplitude and interval time for every region of the brain.
- burst_coincidence_connectivity.m creates a functional connectivity matrix using only coincident bursts between region pairs.
- burst_pdd.m calculates the phase difference derivative over a full second centred on the middle of each burst. This is done for coincident bursts between regions, and also for 'non-coincident bursts' where there is a burst in the seed region, but not the test region.



