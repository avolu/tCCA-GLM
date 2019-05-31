clear all;

% ##### FOLLOWING TWO LINES NEED CHANGE ACCORDING TO USER!
%Alex
path.code = 'D:\Office\Research\Software - Scripts\Matlab\Regression tCCA GLM\tCCA-GLM'; addpath(genpath(path.code)); % code directory
path.dir = 'C:\Users\avolu\Google Drive\tCCA_GLM_PAPER\FB_RESTING_DATA'; % data directory

% #####
%% simulated data file names
filename = 'resting_sim';
%% load ground truth hrf
hrf = load([path.code '\sim HRF\hrf_simdat_100.mat']);

set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig))
set(groot,'defaultAxesCreateFcn',@(ax,~)set(ax.Toolbar,'Visible','off'))
sbjfolder = {'Subj33','Subj34','Subj36','Subj37','Subj38','Subj39', 'Subj40', 'Subj41', 'Subj43', 'Subj44','Subj46','Subj47','Subj49','Subj51'};


%% Options/Parameter Settings
rhoSD_ssThresh = 15;  % mm
flag_save = 0;
flag_conc = 1; % if 1 CCA inputs are in conc, if 0 CCA inputs are in intensity
% results eval parameters
eval_param.HRFmin = -2;
eval_param.HRFmax = 17; % used only for block design runs
eval_param.Hb = 1; % 1 HbO / 0 HbR (for block only)
eval_param.pre = 5;  % HRF range in sec to calculate ttest
eval_param.post = 10;
% CCA parameters
flags.pcaf =  [0 0]; % no pca of X or AUX

%motion artifact detection
motionflag = true;
%plot flag
flag_plot = true;
% flag for mse/corr for each trial (1 = get sum of mse for each trial, 0 = get mse for average estimated hrf)
flag_trial = 0;

% Available parameters
tlags = 0:1:10;
stpsize = 2:2:24;
cthresh = 0:0.1:0.9;

%% Fixed parameters for CCA
pOptfix = [4 8 6];
timelag = pOptfix(1); %loop across timelags
sts = pOptfix(2); % stepsizes
ctr = pOptfix(3); % correlation threshold

tic;

%% Choose subject 
% ALL
% slist = numel(sbjfolder);
% one
slist = 1;

for sbj = 1:slist % loop across subjects
    disp(['subject #' num2str(sbj)]);
       
    % change to subject directory
    cd([path.dir filesep sbjfolder{sbj} filesep]);
    
    %% load data
    [fq, t, AUX, d_long, d_short, d0_long, d0_short, d, d0, SD, s, lstLongAct,lstShortAct,lstHrfAdd] = load_nirs(filename,flag_conc);
    
    %% lowpass filter AUX signals
    AUX = hmrBandpassFilt(AUX, fq, 0, 0.5);
    %% AUX signals
    AUX = [AUX, d0_short]; % full AUX = [acc1 acc2 acc3 PPG BP RESP, d_short];
    %% zscore AUX signals
    AUX = zscore(AUX);
    
    %% check if the number of time points is odd/even, if odd make it even... (number of embedded should be the same)
    if mod(size(AUX,1),2) == 1
        AUX(end,:)=[];
        d(end,:)=[];
        d_long(end,:)=[];
        d_short(end,:)=[];
        d0(end,:)=[];
        d0_long(end,:)=[];
        d0_short(end,:)=[];
        t(end,:)=[];
        s(end,:)=[];
    end
    
    % create data split indices
    len = size(AUX,1);
    spltIDX = {1:len/2,len/2+1:len};
    trntst = {[1,2], [2,1]};
    
    %% run test and train CV splits
    for tt = 1:2
        tstIDX = spltIDX{trntst{tt}(1)};
        trnIDX = spltIDX{trntst{tt}(2)};
        
        %% convert testing fNIRS data to concentration and detect motion artifacts
        dod = hmrIntensity2OD(d(tstIDX,:));
        
        if motionflag
            [tIncAuto] = hmrMotionArtifact(dod,fq,SD,ones(size(d,1),1),0.5,1,30,5);
            [s,tRangeStimReject] = enStimRejection(t(tstIDX,:),s,tIncAuto,ones(size(d,1),1),[-2  10]);
        end
        
        dod = hmrBandpassFilt(dod, fq, 0, 0.5);
        dc{tt} = hmrOD2Conc( dod, SD, [6 6]);
        
        %% Perform GLM with SS
        [yavg_ss, yavgstd_ss, tHRF, nTrialsSS, d_ss, yresid_ss, ysum2_ss, beta_ss, yR_ss] = ...
            hmrDeconvHRF_DriftSS(dc{tt}, s(tstIDX,:), t(tstIDX,:), SD, [], [], [eval_param.HRFmin eval_param.HRFmax], 1, 1, [0.5 0.5], rhoSD_ssThresh, 1, 3, 0);
        
        %% CCA EVAL
        %% set stepsize for CCA
        param.tau = sts; %stepwidth for embedding in samples (tune to sample frequency!)
        param.NumOfEmb = ceil(timelag*fq / sts);
        
        %% Temporal embedding of auxiliary data from testing split
        aux_sigs = AUX(tstIDX,:);
        aux_emb = aux_sigs;
        for i=1:param.NumOfEmb
            aux=circshift( aux_sigs, i*param.tau, 1);
            aux(1:2*i,:)=repmat(aux(2*i+1,:),2*i,1);
            aux_emb=[aux_emb aux];
        end
        
        %% set correlation trheshold for CCA to 0 so we dont lose anything here
        param.ct = 0;   % correlation threshold
        %% Perform CCA on training data % AUX = [acc1 acc2 acc3 PPG BP RESP, d_short];
        % use test data of LD channels without synth HRF
        X = d0_long(trnIDX,:);
        [REG_trn{tt},  ADD_trn{tt}] = perf_temp_emb_cca(X,AUX(trnIDX,:),param,flags);
        disp(['split: ' num2str(tt) ', tlag: ' num2str(timelag) ', stsize: ' num2str(sts) ', ctrhesh: ' num2str(ctr)])
        
        %% now use correlation threshold for CCA outside of function to avoid redundant CCA recalculation
        % overwrite: auxiliary cca components that have
        % correlation > ctr
        compindex=find(ADD_trn{tt}.ccac>ctr);
        %overwrite: reduced mapping matrix Av
        ADD_trn{tt}.Av_red = ADD_trn{tt}.Av(:,compindex);
        
        %% Calculate testig regressors with CCA mapping matrix A from testing
        REG_tst = aux_emb*ADD_trn{tt}.Av_red;
        
        %% Perform GLM with CCA
        [yavg_cca, yavgstd_cca, tHRF, nTrials(tt), d_cca, yresid_cca, ysum2_cca, beta_cca, yR_cca] = ...
            hmrDeconvHRF_DriftSS(dc{tt}, s(tstIDX,:), t(tstIDX,:), SD, REG_tst, [], [eval_param.HRFmin eval_param.HRFmax], 1, 1, [0.5 0.5], 0, 0, 3, 0);
        
        
        %% list of channels with stimulus
        lst_stim = find(s(tstIDX,:)==1);
        if lst_stim(1) < abs(eval_param.HRFmin) * fq
            lst_stim = lst_stim(2:end);
        end
        if size(s(tstIDX,:),1) < lst_stim(end) + abs(eval_param.HRFmax) * fq
            lst_stim = lst_stim(1:end-1);
        end
        
        
        %% EVAL / PLOT
        [DET_SS(:,:,tt), DET_CCA(:,:,tt), pval_SS(:,:,tt), pval_CCA(:,:,tt), ...
            ROCLAB, MSE_SS(:,:,tt), MSE_CCA(:,:,tt), CORR_SS(:,:,tt), CORR_CCA(:,:,tt)] = ...
            results_eval(sbj, d_ss, d_cca, yavg_ss, yavg_cca, tHRF, timelag, sts, ctr, lst_stim, SD, fq, lstHrfAdd, lstLongAct, eval_param, flag_plot, path, hrf, flag_trial, nTrials(tt));
        % Dimensions of output metrics
        % #CH x 2(Hbo+HbR) x 2 (cv split) x tlag x stepsize x corrthres
        % old:  #CH x 2(Hbo+HbR) x 2 (cv split) x SBJ x tlag x stepsize x corrthres
    end
end


toc;



