% Script to calculate the resting state connectivity matrix using only
% coincident bursts.
% The Jaccard index is computed for each pair of regions to create a single
% matrix for each subject. These are then averaged over subjects and
% converted to a pseudo-z-statistic and plotted.

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
for sub = 1:length(sub_IDs)
    clearvars -except sub sub_IDs options jaccard_conn_map
    
    disp(['Subject ',sub_IDs(sub,:),', reading in data'])
    % Virtual electrode timecourses for 78 brain regions.
    % 1-48Hz data:
    VE_dir = strcat('/home/Zelekha/VEs_RSTST/',sub_IDs(sub,:))
    cd(VE_dir)
    load VEf_b_1_48.mat
    % In our case, the VE matrix is (78 brain regions x length of dataset)
    T = size(VEf_b_1_48,2);
    % Downsample from 600 to 100Hz
    [data_new,T_new] = downsampledata(VEf_b_1_48',T,options.Fs,600);
    
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
%% Classify burst state for each region
    for reg = 1:78 % AAL regions
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
        
        % Burst state is the one that correlates best with beta envelope:
        [a, burst_state] = max(corr_tmp);
        
        % Classify bursts by thresholding the state probability timecourse
        % (gamma):
        gamma_thresh = 2/3;
        burst_mask_gamma(reg,:) = Gamma_reg(:,burst_state)>gamma_thresh;
        
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Burst coincidence
    disp(['Subject ',sub_IDs(sub,:),', calculating coincident bursts'])
    for R1 = 1:78
        % Get data for R1
        burst_model_R1 = burst_mask_gamma(R1,:);        
        for R2 = 1:78
            % Get data for R2
            burst_model_R2 = burst_mask_gamma(R2,:);
            % Now count the number of burst overlaps
            jaccard_conn_map(R1,R2,sub) = jaccard(burst_model_R1,burst_model_R2);         
        end
    end   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plots
% No z-score, just jaccard
mean_conn_map = mean(jaccard_conn_map,3); % Average over participants
mean_conn_map(mean_conn_map == 1) = nan; % Set diagonal to nan
figure; imagesc(mean_conn_map); caxis([min(mean_conn_map(:)), 0.27])

% Compute pseudo-z-statistic using mean matrix as baseline
mean_of_mat = nanmean(mean_conn_map(:));
std_of_mat = nanstd(mean_conn_map(:));
z_conn_mat_1 = (mean_conn_map-mean_of_mat)./std_of_mat;
figure; imagesc(z_conn_mat_1); caxis([min(z_conn_mat_1(:)), 4.5])
yticks([5, 14, 25, 37, 44, 53, 64, 76]); 
yticklabels({'L.M.Frontal', 'L.P.Motor', 'L.Calcarine',...
    'L.Cingulum', 'R.M.Frontal', 'R.P.Motor', 'R.Calcarine', 'R.Cingulum'});
xticks([5, 14, 25, 37, 44, 53, 64, 76]); 
xticklabels({'L.M.Frontal', 'L.P.Motor', 'L.Calcarine',...
    'L.Cingulum', 'R.M.Frontal', 'R.P.Motor', 'R.Calcarine', 'R.Cingulum'});
xtickangle(45)
h = colorbar; ylabel(h,'z-score')
set(gca,'Fontsize',16,'LineWidth',2.5,'TickLength',[0.025,0.025]); 



