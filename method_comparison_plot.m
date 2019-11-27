% Script to plot the bursts detected using a simple threshold alongside the
% bursts detected using the time-delay embedded HMM.
% The state probability timecourses and spectra are also plotted. 

% These plots are for a single subject, in a single example region in the
% brain.

% For full details describing how to use the HMM toolbox, go to 
% https://github.com/OHBA-analysis/HMM-MAR/wiki

% Please note that you will have to adapt this script to accomodate your
% local filepaths.

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
%% HMM-AR options
% Set the HMM Options
options = struct(); % Create options struct
no_states = 3;
Hz = 100; % the frequency of the downsampled data
lags = 11; % range between 3 and 11.
options.K = no_states;
options.standardise = 1; % since we hacked the T vector
options.verbose = 1;
options.Fs = Hz;
options.order = 0;
options.embeddedlags = -lags:lags; 
options.zeromean = 1;
options.covtype = 'full'; 
options.useMEX = 1; % runs much faster
options.dropstates = 1;
options.DirichletDiag = 10; % diagonal of prior of trans-prob matrix (default 10 anyway)
options.useParallel = 0;
options.win = 256; % For the multitaper. larger number = smoother spectrum

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load data:
for sub = 1
    clearvars -except sub sub_IDs options
    disp(['Subject ',sub_IDs(sub,:),', reading in data'])
    % Virtual electrode timecourses for 78 brain regions.
    % 1-48Hz data:
    VE_dir = strcat('/home/Zelekha/VEs_RSTST/',sub_IDs(sub,:));
    hmmar_dir = strcat('/home/Zelekha/HMM_RSTST/',sub_IDs(sub,:),'/');  

    % Get 1-48 Hz filtered VE timecourse
    cd(VE_dir)
    load VEf_b_1_48.mat
    T = size(VEf_b_1_48,2);
    % Downsample from 600 to 100Hz
    [data_new,T_new] = downsampledata(VEf_b_1_48',T,100,600);

    % Hack to allow some temporally local variance normalisation
    dummy_trial_length = 4; %secs
    dummy_trial_length = dummy_trial_length*100; % in num of samples
    T_cat = ones(ceil(size(data_new,1)/dummy_trial_length),1)*dummy_trial_length; 
    tmp = rem(size(data_new,1),dummy_trial_length);
    if tmp>0
        T_cat(end)=tmp;
    end
    T_new = T_cat;
    data_new = data_new';

    % Load 3-state HMMAR output
    cd(hmmar_dir)
    load hmm_all_reg.mat % HMM model
    hmm_all_reg_TE = hmm_all_reg;
    load Gamma_all_reg.mat % state probability timecourses
    Gamma_all_reg_TE = Gamma_all_reg;
    
    % Load the 13-33Hz VE data for the simple threshold burst detection
    cd(VE_dir)
    load VEf_b_13_33.mat
    T = size(VEf_b_13_33,2);
    % Downsample from 600 to 100Hz
    [data_new_13_33,T_new_13_33] = downsampledata(VEf_b_13_33',T,options.Fs,600);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Classify burst state
    for reg = 16 % Select from AAL regions
        VE_reg = data_new(reg,:); % single channel data
        data_reg = normalise(VE_reg');

        % Correct Gamma time-course for lags in AR model
        Gamma_reg = Gamma_all_reg_TE{reg};
        Gamma_reg = padGamma(Gamma_reg,T_new,options);
        
        % 13-30Hz data:
        [wt,wf] = cwt(data_reg,'amor',100); % morlet wavelet
        % NOTE. Need to be using at least matlab v 2016a for the cwt
        % function.
        beta_wavelet_freqs = wf > 13 & wf < 30;
        % Envelope of beta oscillations
        HVEf_b_reg = mean(abs(wt(beta_wavelet_freqs,:)),1)';
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Correlation of burst probability timecourse with beta envelope
        corr_tmp = [];
        for k = 1:size(Gamma_reg,2)
            corr_tmp(k) = corr(Gamma_reg(:,k),HVEf_b_reg);
        end
        
        % The burst state is the one that correlates most with the beta
        % envelope
        [a burst_state] = max(corr_tmp);

        % Classify burst segments by thresholding the probability
        % timecourse (gamma):
        gamma_thresh = 2/3;
        burst_mask_gamma(reg,:) = Gamma_reg(:,burst_state)>gamma_thresh;
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% State spectra
        T = length(data_reg);
        fit = hmmspectramt(data_reg,T,Gamma_reg,options); % statewise multitaper
        figure; plot(fit.state(2).f,fit.state(2).psd,'-','LineWidth',2);
        hold on; plot(fit.state(1).f,fit.state(1).psd,'-','LineWidth',2);
        hold on; plot(fit.state(3).f,fit.state(3).psd,'-','LineWidth',2);
        xlabel('Frequency (Hz)'); ylabel('P.S.D.');
        title(['Burst state: ',num2str(burst_state)]);
        set(gca,'FontSize',16); legend('state 1', 'state 2', 'state 3');
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Comparison    
        %threshold = 2.5*std:
        data_reg_13_33 = data_new_13_33(:,reg);
        % Envelope of beta oscillations
        HVEf_b = abs(hilbert(data_reg_13_33));
        std_reg = std(data_reg_13_33);
        bursts_SD25 = (HVEf_b > 2.5*std_reg);
        burst_vals_SD25 = HVEf_b.*bursts_SD25;
        burst_vals_SD25(burst_vals_SD25 == 0) = nan;
        
        % HMM-TE burst categorisation:
        burst_vals_HMM_AR = data_reg'.*burst_mask_gamma(reg,:);
        burst_vals_HMM_AR(burst_vals_HMM_AR == 0) = nan;
        
        start_time = 20*options.Fs; % 20s
        end_time = 40*options.Fs; % 40s
        time_s = (start_time:end_time)/100;
        figure;
        % Threshold
        subplot(4,1,1); plot(time_s,HVEf_b(start_time:end_time)); grid on; hold on;
        plot(time_s,burst_vals_SD25(start_time:end_time));
        str = "Hilbert" + newline + "Amplitude" + newline + "(nAm)";
        xlabel('Time (s)'); ylabel(str); set(gca,'Fontsize',20)
        % HMM-TE
        subplot(4,1,2); plot(time_s,data_reg(start_time:end_time)); grid on; hold on;
        plot(time_s,burst_vals_HMM_AR(start_time:end_time));
        str = "Regional" + newline + "Timecourse" + newline + "(A.U.)";
        xlabel('Time (s)'); ylabel(str); set(gca,'Fontsize',20)
        % Gamma timecourse
        subplot(4,1,3); area(time_s,Gamma_reg(start_time:end_time,burst_state));
        grid on; hold on; 
        plot(time_s,burst_mask_gamma(reg,start_time:end_time),'-','LineWidth',2);
        ylim([0, 1.1]); str = "State" + newline + "Probability";
        xlabel('Time (s)'); ylabel(str); set(gca,'Fontsize',20)
        
        % Wavelet
        subplot(4,1,4); tbins = start_time:end_time;
        contourf( tbins/100,wf, sqrt(abs(wt(:,tbins))), 36 ,'linestyle','none' )
        grid on;hold on
        str = "Frequency";
        str = str + newline + "(Hz)";
        xlabel('Time (s)'); ylabel(str);   
        set(gca,'Fontsize',20)
    end
end



