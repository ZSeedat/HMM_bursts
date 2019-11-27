% Script to calculate the burst features for the HMM_TE bursts.
% Features include number of bursts, burst power, burst duration, burst 
% state correlation and time between bursts.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
clear all;
% File with the hmm toolbox
addpath('/home/Zelekha/')
addpath(genpath('/home/Zelekha/HMM_MAR_Master/HMM-MAR-master/'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Filenames for UKMP data - resting state
load sub_IDs.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HMM Options:
options = struct(); % Create options struct
no_states = 3;
Hz = 100; % the frequency of the downsampled data
options.K = no_states;
options.standardise = 1; % since we hacked the T vector
options.verbose = 1;
options.Fs = Hz;
options.order = 0;
lags = 11; % a sensible range is between 3 and 11.
options.embeddedlags = -lags:lags; 
options.zeromean = 1;
options.covtype = 'full'; 
options.useMEX = 1; % runs much faster
options.dropstates = 1;
options.DirichletDiag = 10; % diagonal of prior of trans-prob matrix (default 10 anyway)
options.useParallel = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load data:
for sub = 1:length(sub_IDs);
    clearvars -except sub sub_IDs options num_bursts burst_pow burst_dur corr_vals burst_ISL burst_overlap_all
    
    disp(['Subject ',sub_IDs(sub,:),', reading in data'])
    % Virtual electrode timecourses for 78 brain regions.
    % 1-48Hz data:
    VE_dir = strcat('/home/Zelekha/VEs_RSTST/',sub_IDs(sub,:));
    cd(VE_dir)
    load VEf_b_1_48.mat
    T = size(VEf_b_1_48,2);
    % Downsample from 600 to 100Hz
    [data_new,T_new] = downsampledata(VEf_b_1_48',T,100,600);

    % Hack to allow some temporally local variance normalisation
    dummy_trial_length = 4; % secs
    dummy_trial_length = dummy_trial_length*100; % in num of samples
    T_cat = ones(ceil(size(data_new,1)/dummy_trial_length),1)*dummy_trial_length; 
    tmp = rem(size(data_new,1),dummy_trial_length);
    if tmp > 0
        T_cat(end) = tmp;
    end
    T_new = T_cat;
    data_new = data_new';

    % Load the HMM-TE output
    HMM_dir = strcat('/home/Zelekha/HMM_RSTST/',sub_IDs(sub,:),'/');
    cd(HMM_dir)
    load hmm_all_reg.mat
    hmm_all_reg_TE = hmm_all_reg;
    load Gamma_all_reg.mat 
    Gamma_all_reg_TE = Gamma_all_reg;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Classify burst state
    for reg = 1:78 % 78 AAL regions 
        VE_reg = data_new(reg,:); % single channel data
        data_reg = normalise(VE_reg');

        % Correct Gamma time-course for lags in AR model
        Gamma_reg = Gamma_all_reg_TE{reg};
        Gamma_reg = padGamma(Gamma_reg,T_new,options);
        
        % 13-30Hz data:
        [wt,wf] = cwt(data_reg,'amor',100); % morlet wavelet
        beta_wavelet_freqs = wf > 13 & wf < 30;
        % Envelope of beta oscillations
        HVEf_b_reg = mean(abs(wt(beta_wavelet_freqs,:)),1)';
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Correlation of burst probability timecourse with beta envelope
        corr_tmp = [];
        for k = 1:size(Gamma_reg,2)
            corr_tmp(k) = corr(Gamma_reg(:,k),HVEf_b_reg);
        end

        [a burst_state] = max(corr_tmp);
        corr_vals(sub,reg) = max(corr_tmp);

        % Classify bursts using thresholded gamma:
        gamma_thresh = 2/3;
        burst_mask_gamma(reg,:) = Gamma_reg(:,burst_state)>gamma_thresh;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Burst duration
        % Burst lifetimes
        LTs_hmmar = getStateLifeTimes(Gamma_reg,size(Gamma_reg,1),options,0,gamma_thresh);
        burst_LTs = LTs_hmmar{burst_state};
        burst_dur(sub,reg) = mean(burst_LTs)*1000/options.Fs; % milliseconds

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%% Number of bursts
        num_bursts(sub,reg) = length(burst_LTs)/(size(Gamma_reg,1)/options.Fs);
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Burst amplitude
        burst_envelope = HVEf_b_reg'.*burst_mask_gamma(reg,:);
        burst_starts = find(diff(burst_mask_gamma(reg,:)) == 1);
        burst_ends = find(diff(burst_mask_gamma(reg,:)) == -1);
        if burst_starts(1) > burst_ends(1)
            burst_ends(1) = [];
        end
        if burst_starts(end) > burst_ends(end)
            burst_starts(end) = [];
        end
        clear burst_max
        for b = 1:length(burst_starts)
            burst_max(b) = max(burst_envelope(burst_starts(b):burst_ends(b)));
        end
        burst_pow(sub,reg) = mean(burst_max);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Time between bursts
        state_ITs = getStateIntervalTimes(Gamma_reg,size(Gamma_reg,1),options,0,gamma_thresh);
        burst_ITs = state_ITs{burst_state};
        burst_ISL(sub,reg) = mean(burst_ITs)*1000/options.Fs; % milliseconds
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save output
save_to = '/home/Zelekha/HMM_RSTST/';
cd(save_to)
save burst_dur.mat burst_dur
save num_bursts.mat num_bursts
save burst_pow.mat burst_pow
save burst_ISL.mat burst_ISL


