% Script to investigate the coherence of coincident bursts in two 
% strongly-connected regions. Uses phase difference derivative (PDD).
% PDD code is from Lucrezia Liuzzi.

% The coincident bursts and the non-coincident bursts (where there is a
% burst in the seed region, but not the test region) are compared in a plot
% of PDD at the end of the script.

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
pdd_all_c = []; % coincident bursts
pdd_all_nc = []; % non-coincident bursts
num_c = 0; % count the c bursts
num_nc = 0; % count the nc bursts
for sub = 1:length(sub_IDs);
    clearvars -except sub sub_IDs options pdd_all_c pdd_all_nc num_c num_nc
    
    disp(['Subject ',sub_IDs(sub,:),', reading in data'])
    % Virtual electrode timecourses for 78 brain regions.
    % 1-48Hz data:
    VE_dir = strcat('/home/Zelekha/VEs_RSTST/',sub_IDs(sub,:));
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
    load hmm_all_reg_TE.mat
    hmm_all_reg_TE = hmm_all_reg;
    load Gamma_all_reg_TE.mat
    Gamma_all_reg_TE = Gamma_all_reg;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Classify burst state
    % Relevent regions only
    Reg_1 = [14,16,17,22,10,3]; % seed region (LMotor,LSensory,LParietal,LVisual,LFrontal,LF_orb)
    Reg_2 = [53,55,56,61,49,42]; % connected region (RM,RS,RP,RV,RF,RF_orb)
    for connection = 6 % ONLY RUN ONE CONNECTION AT A TIME OR IT WILL OVER-WRITE IT
        for reg = [Reg_1,Reg_2] % AAL regions
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
            
            % Select burst state as the one that correlates best with the
            % beta envelope
            [a burst_state] = max(corr_tmp);
            
            % Classify bursts by thresholding the state probability
            % timecourse (gamma)
            gamma_thresh = 2/3;
            burst_mask_gamma(reg,:) = Gamma_reg(:,burst_state)>gamma_thresh;            
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Burst coherence
        % Look at the phase +/- 0.5s either side of a burst centre.
        % Calculate the phase coherence for coincident bursts vs non-coincident
        % bursts.
        disp(['Subject ',sub_IDs(sub,:),', calculating pdd'])
        R1 = Reg_1(connection); % Seed region
        clear burst_starts burst_ends
        % Get data for R1
        VE_R1 = data_new(R1,:);
        data_R1 = normalise(VE_R1');
        % Calculate phase
        phi_R1 = unwrap(angle(data_R1));
        burst_model_R1 = burst_mask_gamma(R1,:);
        % Find bursts
        burst_starts = find(diff(burst_model_R1) == 1);
        burst_ends = find(diff(burst_model_R1) == -1);
        
        % Align the burst starts and ends
        if burst_starts(1) > burst_ends(1)
            burst_ends(1) = [];
        end
        if burst_starts(end) > burst_ends(end)
            burst_starts(end) = [];
        end
        
        % Find the index of each burst centre
        burst_mids_R1 = burst_starts + (burst_ends - burst_starts)/2;
        
        R2 = Reg_2(connection); % Connected region
        clear burst_starts burst_ends
        % Get data for R2
        VE_R2 = data_new(R2,:);
        data_R2 = normalise(VE_R2');
        % Calculate phase
        phi_R2 = unwrap(angle(data_R2));
        burst_model_R2 = burst_mask_gamma(R2,:);
        % Find bursts
        burst_starts = find(diff(burst_model_R2) == 1);
        burst_ends = find(diff(burst_model_R2) == -1);
        
        % Align the burst starts and ends
        if burst_starts(1) > burst_ends(1)
            burst_ends(1) = [];
        end
        if burst_starts(end) > burst_ends(end)
            burst_starts(end) = [];
        end
        
        % Find the index of each burst centre in region 2
        burst_mids_R2 = burst_starts + (burst_ends - burst_starts)/2;
        
        % Now find the centres of the burst overlaps
        coinc_bursts = burst_model_R1.*burst_model_R2;
        overlap_starts = find(diff(coinc_bursts) == 1);
        overlap_ends = find(diff(coinc_bursts) == -1);
        
        % Align the burst starts and ends
        if overlap_starts(1) > overlap_ends(1)
            overlap_ends(1) = [];
        end
        if overlap_starts(end) > overlap_ends(end)
            overlap_starts(end) = [];
        end
        
        % Find the centre of each overlapping section
        overlap_mids = overlap_starts + (overlap_ends - overlap_starts)/2;
        % Remove any indices closer than 0.5s to the start/end of
        % data for the overlapping bursts
        overlap_mids(overlap_mids < 0.5*options.Fs+1) = [];
        overlap_mids(overlap_mids > length(data_R1)-0.5*options.Fs) = [];
        overlap_mids = round(overlap_mids);
        
        % Remove any indices closer than 0.5s to the start/end of
        % data for any of the bursts in R1
        burst_mids_R1(burst_mids_R1 < 0.5*options.Fs+1) = [];
        burst_mids_R1(burst_mids_R1 > length(data_R1)-0.5*options.Fs) = [];
        burst_mids_R1 = round(burst_mids_R1);
        % ...similarly for R2
        burst_mids_R2(burst_mids_R2 < 0.5*options.Fs+1) = [];
        burst_mids_R2(burst_mids_R2 > length(data_R2)-0.5*options.Fs) = [];
        burst_mids_R2 = round(burst_mids_R2);
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Coincident bursts
        % Find the bursts resposible for the overlaps in both regions
        for b = 1:length(overlap_mids) % for each overlapping burst
            % Region 1:
            min_R1_diff = min(abs(burst_mids_R1-overlap_mids(b)));
            if min_R1_diff == 0 % if the burst centres overlap exactly
                R1_ind = find(abs(burst_mids_R1-overlap_mids(b)) == 0);
            else
                R1_ind = find(abs(burst_mids_R1-overlap_mids(b)) == min_R1_diff);
            end
            R1_ind = R1_ind(1); % where bursts are equidistant from the overlap centre, take the first one
            R1_mid(b) = burst_mids_R1(R1_ind);
            
            % Now chop out each burst responsible for an overlap
            phi_segments_R1_c(b,:) = phi_R1(R1_mid(b)-0.5*options.Fs:R1_mid(b)+0.5*options.Fs);
            
            % Region 2:
            min_R2_diff = min(abs(burst_mids_R2-overlap_mids(b)));
            if min_R2_diff == 0 % if the burst centres overlap exactly
                R2_ind = find(abs(burst_mids_R2-overlap_mids(b)) == 0);
            else
                R2_ind = find(abs(burst_mids_R2-overlap_mids(b)) == min_R2_diff);
            end
            R2_ind = R2_ind(1); % where bursts are equidistant from the overlap centre, take the first one only
            R2_mid(b) = burst_mids_R2(R2_ind);
            
            % Now chop out each burst responsible for an overlap
            phi_segments_R2_c(b,:) = phi_R2(R2_mid(b)-0.5*options.Fs:R2_mid(b)+0.5*options.Fs);
        end
        
        % Number coincident bursts:
        num_c_sub = b;
        num_c = num_c + b;
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Non-coincident bursts
        % Remove coincident burst centres from R1
        burst_mids_R1_nc = burst_mids_R1;
        for b = 1:length(R1_mid)
            burst_mids_R1_nc(burst_mids_R1_nc == R1_mid(b)) = [];
        end
        for b = 1:length(burst_mids_R1_nc) % for the bursts left
            % Region 1:
            mid_b = burst_mids_R1_nc(b);
            
            % Now chop out each non-coincident burst from the seed region
            phi_segments_R1_nc(b,:) = phi_R1(mid_b-0.5*options.Fs:mid_b+0.5*options.Fs);
            
            % Select analogous data from connected region
            phi_segments_R2_nc(b,:) = phi_R2(mid_b-0.5*options.Fs:mid_b+0.5*options.Fs);
        end
        
        % Number non-coincident bursts:
        num_nc_sub = b;
        num_nc = num_nc + b;
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Phase Difference Derivative
        % Make sure we match the number of coincident and non-coincident
        % bursts
        if num_c_sub > num_nc_sub
            num_diff = num_c_sub - num_nc_sub;
            phi_segments_R1_c(end-num_diff+1:end,:) = [];
            phi_segments_R2_c(end-num_diff+1:end,:) = [];
        else
            num_diff = num_nc_sub - num_c_sub;
            phi_segments_R1_nc(end-num_diff+1:end,:) = [];
            phi_segments_R2_nc(end-num_diff+1:end,:) = [];
        end
            
        % Coincident bursts
        % Phase difference (units [radians])
        PDc = phi_segments_R1_c-phi_segments_R2_c;
        % Derivative (units [radians/sample])
        PDDc = diff(PDc,1,2); % Along second dimension
        % Convert into exponential function
        pdd_exp_c = exp(-abs(PDDc));  % units [ radians/sample]
        % Average so that we have a single pdd value per subject
        mean_pdd_c_sub = mean(pdd_exp_c,1);
        
        % Non-coincident bursts
        % Phase difference (units [radians])
        PDnc = phi_segments_R1_nc-phi_segments_R2_nc;
        % Derivative (units [radians/sample])
        PDDnc = diff(PDnc,1,2); % Along second dimension
        % Convert into exponential function
        pdd_exp_nc = exp(-abs(PDDnc));  % units [ radians/sample]
        % Average so that we have a single pdd value per subject
        mean_pdd_nc_sub = mean(pdd_exp_nc,1);
    end
    
    % Concatenate average pdd value for each subject
    pdd_all_c = cat(1,pdd_all_c,mean_pdd_c_sub);
    % Concatenate non-coincident pdd for each subject
    pdd_all_nc = cat(1,pdd_all_nc,mean_pdd_nc_sub);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot
% Coincident bursts
mean_pdd = mean(pdd_all_c,1);
% Smooth function a little
smooth_pdd = smooth(mean_pdd,5);
% Set error boundaries
err_pdd = std(pdd_all_c,1)/sqrt(size(pdd_all_c,1));
x = 1:length(smooth_pdd);
lo = smooth_pdd' - err_pdd;
hi = smooth_pdd' + err_pdd;
figure;
% Define patch
hp1 = patch([x, x(end:-1:1)], [lo, hi(end:-1:1)], 'r');
hold on;
% Define line
hl = plot(smooth_pdd);
set(hp1, 'facecolor', [173 216 230]/255, 'edgecolor', 'none');
set(hl, 'color', [65 105 225]/255);
title('L-R Oribitofrontal Cortex')
xticks([10.*(0:10)]); xticklabels({'-0.5','-0.4','-0.3','-0.2','-0.1','0','0.1','0.2','0.3','0.4','0.5'});
xlabel('time, s'); ylabel('PDD'); set(gca,'fontsize',16);

% Non-coincident bursts
mean_pdd_nc = mean(pdd_all_nc,1);
% Smooth function a little
smooth_pdd_nc = smooth(mean_pdd_nc,5);
% Set error boundaries
err_pdd_nc = std(pdd_all_nc,1)/sqrt(size(pdd_all_nc,1));
lo_nc = smooth_pdd_nc' - err_pdd_nc;
hi_nc = smooth_pdd_nc' + err_pdd_nc;
% Define patch
hp2 = patch([x, x(end:-1:1)], [lo_nc, hi_nc(end:-1:1)], 'r');
hold on;
% Define line
h2 = plot(smooth_pdd_nc);
set(hp2, 'facecolor', [255 218 185]/255, 'edgecolor', 'none');
set(h2, 'color', [255 69 0]/255);
alpha(0.4); % transparency value

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate p-values
% Between -0.2 and 0.2 of a second
for_p_c = mean(pdd_all_c(:,30:70),2); % coincident bursts
for_p_nc = mean(pdd_all_nc(:,30:70),2); % non-coincident bursts
p = ranksum(for_p_c,for_p_nc);






