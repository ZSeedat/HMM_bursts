% Script to load the 1-48 Hz VE timecourse for each region, followed
% by running a three-state HMM with temporal embedding (HMM-TE) on each 
% region separately. 
% HMM_TE output is then saved for future analysis.

% For full details describing how to use the HMM toolbox, go to 
% https://github.com/OHBA-analysis/HMM-MAR/wiki

% Please note that you will have to adapt this script to accomodate your
% local filepaths.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
clear all; close all;
% Add filepath with the HMM_Master file in it
addpath('/home/Zelekha/')
addpath(genpath('/home/Zelekha/HMM_MAR_Master/HMM-MAR-master/'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Filenames for UKMP data - resting state
load sub_IDs.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run HMM_TE on each subject and region individually
for sub = 1:length(sub_IDs)
    clear global
    display(['Reading in data from subject ',sub_IDs(sub,:)]);
    
    % Get 1-48 Hz filtered VE timecourse
    cd(strcat('/home/Zelekha/VEs_RSTST/',sub_IDs(sub,:)));
    load VEf_b_1_48.mat
    % In our case, the VE matrix is (78 brain regions x length of dataset)
    T = length(VEf_b_1_48);
    samp_freq = 600; % sampling frequency, Hz
    % Downsample from 600 to 100Hz
    new_freq = 100;
    [data_new,T_new] = downsampledata(VEf_b_1_48',T,new_freq,samp_freq);
    
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
    
    % Set the HMM Options
    Hz = 100; % the frequency of the downsampled data
    lags = 11; % sensible to range between 3 and 11.
    no_states = 3;
    options = struct(); % Create options struct
    options.K = no_states;
    options.standardise = 1;
    options.verbose = 1;
    options.Fs = Hz;
    options.order = 0;
    options.embeddedlags = -lags:lags; 
    options.zeromean = 1;
    options.covtype = 'full'; 
    options.useMEX = 1; % runs much faster with the compiled mex files
    options.dropstates = 1; % the HMM can drop states if there isn't sufficient evidence to justify this number of states.
    options.DirichletDiag = 10; % diagonal of prior of trans-prob matrix (default 10 anyway)
    options.useParallel = 0;

    % HMM computation, one region at a time
    for reg = 1:78 % 78 AAL Atlas locations
        disp(['Calculating HMM output for subject ',num2str(sub), ', region ',num2str(reg)])
        VEf_b_reg = data_new(reg,:); % single channel data
        data_reg = normalise(VEf_b_reg'); % normalise!!
        [hmm, Gamma] = hmmmar(data_reg,T_new,options); % hmm inference
        hmm_all_reg{reg} = hmm;
        Gamma_all_reg{reg} = Gamma;
    end
    
    % Save the output
    save_to = strcat('/home/Zelekha/HMM_RSTST/',sub_IDs(sub,:),'/');
    if exist(save_to);                                                                                       
        cd(save_to); 
    else
        mkdir(save_to);                                                                                      
        cd(save_to);                                                                                         
    end
    save hmm_all_reg.mat
    save Gamma_all_reg.mat
end

