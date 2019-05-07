clear all;

malexflag = 1;
if malexflag
    %Meryem
    path.code = 'C:\Users\mayucel\Documents\PROJECTS\CODES\tCCA-GLM'; addpath(genpath(path.code)); % code directory
    path.dir = 'C:\Users\mayucel\Google Drive\tCCA_GLM_PAPER\FB_RESTING_DATA'; % data directory
    path.save = 'C:\Users\mayucel\Google Drive\tCCA_GLM_PAPER'; % save directory
else
    %Alex
    path.code = 'D:\Office\Research\Software - Scripts\Matlab\Regression tCCA GLM\tCCA-GLM'; addpath(genpath(path.code)); % code directory
    path.dir = 'C:\Users\avolu\Google Drive\tCCA_GLM_PAPER\FB_RESTING_DATA'; % data directory
    path.save = 'C:\Users\avolu\Google Drive\tCCA_GLM_PAPER'; % save directory
end

% #####
filename = 'resting.nirs';
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

% prune
flag_prune = false;
%plot flag
flag_plot = false;


% Optimum CCA parameters (obtained through cross-validation analysis of resting data -see tCCA folder)
timelag = 3;
sts = 16;
cthresh = 0.5;

tic;


for sbj = 1:numel(sbjfolder) % loop across subjects
    disp(['subject #' num2str(sbj)]);
    
    
    %% load data
    load([path.dir filesep sbjfolder{sbj} filesep filename],'-mat');
    
    % sampling frequency
    fq = abs(1/(t(1)-t(2)));
    
    % auxilliary channels
    AUX = aux(:,2:7);
    % acc1 = aux(:,2); % accelerometer
    % acc2 = aux(:,3); % accelerometer
    % acc3 = aux(:,4); % accelerometer
    % PPG = aux(:,5); % pulse
    % BP = aux(:,6); % blood pressure waveform
    % RESP = aux(:,7); % respiration
    
    % lpf data entering CCA
    foo = hmrBandpassFilt(d, fq, 0, 0.5);
    
    
    %% Prune noisy channels
    if flag_prune
        SD = enPruneChannels(d,SD,ones(size(d,1),1),[10000  10000000],5,[0  45],0);
    end
    
    %% Definition channel groups
    ml = SD.MeasList;
    lst = find(ml(:,4)==1); % 690
    rhoSD = zeros(length(lst),1);
    posM = zeros(length(lst),3);
    for iML = 1:length(lst)
        rhoSD(iML) = sum((SD.SrcPos(ml(lst(iML),1),:) - SD.DetPos(ml(lst(iML),2),:)).^2).^0.5;
        posM(iML,:) = (SD.SrcPos(ml(lst(iML),1),:) + SD.DetPos(ml(lst(iML),2),:)) / 2;
    end
    lstLongAct = lst(find(rhoSD> rhoSD_ssThresh & SD.MeasListAct(lst)==1)); % list of long and active channels
    lstShortAct = lst(find(rhoSD< rhoSD_ssThresh & SD.MeasListAct(lst)==1)); % list of long and active channels
    % get d_long and short(and active)
    d_long = [foo(:,lstLongAct), foo(:,lstLongAct+size(d,2)/2)]; % first half 690 nm; second half 830 nm
    d_short = [foo(:,lstShortAct), foo(:,lstShortAct+size(d,2)/2)]; % first half 690 nm; second half 830 nm
    
    
    %% lowpass filter AUX signals
    AUX = hmrBandpassFilt(AUX, fq, 0, 0.5);
    %% AUX signals
    AUX = [AUX, d_short]; % full AUX = [acc1 acc2 acc3 PPG BP RESP, d_short];
    %% zscore AUX signals
    AUX = zscore(AUX);
    
    
    
    
    %% Perform GLM with SS
    [yavg_ss, yavgstd_ss, tHRF, nTrialsSS, d_ss, yresid_ss, ysum2] = hmrDeconvTB_SS3rd(d, s, t, SD, [], [eval_param.HRFmin eval_param.HRFmax],0.5, 0.5, rhoSD_ssThresh);
    
    d_ss(:,[lstLong; lstLong+size(d,2)/2]) = d_ss(:,[lstLong; lstLong+size(d,2)/2])  +  repmat(mean_d_long, size(d_long,1),1);
    
    
    
    
    %% CCA EVAL
    param.tau = sts; %stepwidth for embedding in samples (tune to sample frequency!)
    param.NumOfEmb = ceil(timelag*fq / sts);
    param.ct = cthresh;   % correlation threshold
    
    mean_d_long = mean(d_long,1);
    X = d_long - repmat(mean_d_long, size(d_long,1),1);
    
    
    [REG,  ADD] = perf_temp_emb_cca(X,AUX,param,flags);
    
    
    %% Perform GLM with CCA
    [yavg_cca, yavgstd_cca, tHRF, nTrials, d_cca, yresid_cca, ysum2] = hmrDeconvTB_SS3rd(d, s, t, SD, REG, [eval_param.HRFmin eval_param.HRFmax],0.5, 0.5, 0);
    
    d_cca(:,[lstLong; lstLong+size(d,2)/2]) = d_cca(:,[lstLong; lstLong+size(d,2)/2])  +  repmat(mean_d_long, size(d_long,1),1);
    
    
    x = 1;
    %% EVAL / PLOT
    %     [pval_SS(:,:,sbj), pval_CCA(:,:,sbj)] = ...
    %         results_eval_block(sbj, d_ss, d_cca, yavg_ss, yavg_cca, tHRF, timelag, sts, cthresh, lst_stim, SD, fq, eval_param, flag_plot, path);
    % Dimensions of output metrics
    %
    
end


% 
% %% save data for subject
% disp(['saving...'])
% save([path.save '\Visual_Paradigm_results' '\p_val_results.mat'], 'pval_SS', 'pval_CCA');

