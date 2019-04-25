clear all

%% +++++++++++++++++++++++
%% SCRIPT CONFIGURATION
% +++++++++++++++++++++++
% user: 1 Meryem | 0 Alex
melexflag = 0;
% select which hrf amplitude data: 50 or 100
hrfamp = 100;
% Use only true positives for evaluation of metrics
TP_flag = true;
% number of contours in contour plots
cntno = 10;
% use mean (1) or median (2) in metric contour plots
mflag = 2;
% plot pvalue results or not
pvalflag = false;
% plot other metrics
plotmetrics = false;
%% parameters for determining optima
% normalize metrics: 1 X/max | 2 (X-min)/(max-min)
Jparam.nflag = 2;
% smoothing / optimization metrics: 1 mean, 2 median or 3 all channels
Jparam.mtype = 3;
% Objective function J weights
Jparam.fact.corr = 1;
Jparam.fact.mse =2;
Jparam.fact.pval =0;
Jparam.fact.fscore=2;
Jparam.fact.HbO=1;
Jparam.fact.HbR=1;
% segmentation approach: threshold for segmentation
Jparam.thresh = 0.5;

%% Data
% ##### FOLLOWING TWO LINES NEED CHANGE ACCORDING TO USER!
if melexflag
    %Meryem
    path.code = 'C:\Users\mayucel\Documents\PROJECTS\CODES\tCCA-GLM'; addpath(genpath(path.code)); % code directory
    path.dir = 'C:\Users\mayucel\Google Drive\tCCA_GLM_PAPER\FB_RESTING_DATA'; % data directory
    path.save = 'C:\Users\mayucel\Google Drive\tCCA_GLM_PAPER'; % save directory
    path.cvres50 = 'C:\Users\mayucel\Google Drive\tCCA_GLM_PAPER\CV_results_data_50'; % save directory
    path.cvres100 = 'C:\Users\mayucel\Google Drive\tCCA_GLM_PAPER\CV_results_data_100'; % save directory
else
    %Alex
    path.code = 'D:\Office\Research\Software - Scripts\Matlab\Regression tCCA GLM\tCCA-GLM'; addpath(genpath(path.code)); % code directory
    path.dir = 'C:\Users\avolu\Google Drive\tCCA_GLM_PAPER\FB_RESTING_DATA'; % data directory
    path.save = 'C:\Users\avolu\Google Drive\tCCA_GLM_PAPER'; % save directory
    path.cvres50 = 'C:\Users\avolu\Google Drive\tCCA_GLM_PAPER\CV_results_data_50'; % save directory
    path.cvres100 = 'C:\Users\avolu\Google Drive\tCCA_GLM_PAPER\CV_results_data_100'; % save directory
end

% #####
filename = 'resting_sim';
set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig))
set(groot,'defaultAxesCreateFcn',@(ax,~)set(ax.Toolbar,'Visible','off'))
sbjfolder = {'Subj33','Subj34','Subj36','Subj37','Subj38','Subj39', 'Subj40', 'Subj41', 'Subj43', 'Subj44','Subj46','Subj47','Subj49','Subj51'};

% Validation parameters
tlags = 0:1:10;
stpsize = 2:2:24;
cthresh = 0:0.1:0.9;
evparams.tlags = tlags;
evparams.stpsize = stpsize;
evparams.cthresh = cthresh;

%% load results data from all subjects
% Dimensions of output metrics
% # of sbjs x #CH x 2(Hbo+HbR) x 2 (cv split) x tlag x stepsize x corrthres
for sbj = 1:numel(sbjfolder)
    switch hrfamp
        case 50
            res{sbj} = load([path.cvres50 '\results_sbj' num2str(sbj) '.mat']);
        case 100
            res{sbj} = load([path.cvres100 '\results_sbj' num2str(sbj) '.mat']);
    end
    
    %% append subject matrices here
    CORR_CCA(sbj,:,:,:,:,:,:,:) = res{sbj}.CORR_CCA;
    CORR_SS(sbj,:,:,:,:,:,:,:) = res{sbj}.CORR_SS;
    DET_CCA(sbj,:,:,:,:,:,:,:) = res{sbj}.DET_CCA;
    DET_SS(sbj,:,:,:,:,:,:,:) = res{sbj}.DET_SS;
    MSE_CCA(sbj,:,:,:,:,:,:,:) = res{sbj}.MSE_CCA;
    MSE_SS(sbj,:,:,:,:,:,:,:) = res{sbj}.MSE_SS;
    pval_CCA(sbj,:,:,:,:,:,:,:) = res{sbj}.pval_CCA;
    pval_SS(sbj,:,:,:,:,:,:,:) = res{sbj}.pval_SS;
    nTrials(sbj,:,:,:,:) = res{sbj}.nTrials;
    
end

% true positive only flag
if TP_flag
    pval_SS(find(DET_SS ~= 1)) = NaN;
    pval_CCA(find(DET_CCA ~= 1)) = NaN;
end


%% # of TP/FP/FN/TN channels
% for SS and CCA methods:
foo_SS = permute(DET_SS,[2 1 3 4 5 6 7]);
foo_SS = reshape(foo_SS, size(foo_SS,1), size(foo_SS,2)*size(foo_SS,3)*size(foo_SS,4)*size(foo_SS,5)*size(foo_SS,6)*size(foo_SS,7));
foo_CCA = permute(DET_CCA,[2 1 3 4 5 6 7]);
foo_CCA = reshape(foo_CCA, size(foo_CCA,1), size(foo_CCA,2)*size(foo_CCA,3)*size(foo_CCA,4)*size(foo_CCA,5)*size(foo_CCA,6)*size(foo_CCA,7));
% ROCLAB.name = {'TP','FP','FN','TN', 'PRND'};
for i = 1:size(foo_SS,2)
    % SS
    Ch_TP_SS(i) = sum(foo_SS(:,i)==1);
    Ch_FP_SS(i) = sum(foo_SS(:,i)==-1);
    Ch_FN_SS(i) = sum(foo_SS(:,i)==2);
    Ch_TN_SS(i) = sum(foo_SS(:,i)==-2);
    % CCA
    Ch_TP_CCA(i) = sum(foo_CCA(:,i)==1);
    Ch_FP_CCA(i) = sum(foo_CCA(:,i)==-1);
    Ch_FN_CCA(i) = sum(foo_CCA(:,i)==2);
    Ch_TN_CCA(i) = sum(foo_CCA(:,i)==-2);
end

%% Calculate True/false positive/negative rates, precision, recall, ...
tf_errors


%% Find Global topology and optimum with objective function, includes segmentation approach
% calculate objective function output for all input tupel
fval = J_opt(CORR_CCA, MSE_CCA, pval_CCA, F_score_CCA, Jparam);
% find optimal parameter set
[t,s,c] = ind2sub(size(fval),find(fval == min(fval(:))));
pOpt = [t s c];
disp('=================================================================')
disp(['these parameters minimize the objective function: timelag: ' ...
    num2str(tlags(t)) 's, stepsize: ' num2str(stpsize(s)) 'smpl, corr threshold: ' num2str(cthresh(c))] )
disp('=================================================================')

%% Calculate median/mean metrics
[CORR_CCA,MSE_CCA,pval_CCA,F_score_CCA] = medmean(CORR_CCA, MSE_CCA, pval_CCA, F_score_CCA, mflag);

%% create combined surface plots (depict objective function)
hblab = {'HbO', 'HbR'};

%% Plot Objective function results
% normalize fval
fval = (fval-min(fval(:)))/(max(fval(:))-min(fval(:)));
ttl= 'Objective Function';
contour_plots(fval, ttl, evparams, pOpt, cntno, 'min');

if plotmetrics
    %% plot correlation
    %HBO and HbR
    for hh = 1:2
        ttl= [hblab{hh} ' Correlation'];
        contour_plots(squeeze(CORR_CCA(:,hh,:,:,:)), ttl,evparams, pOpt, cntno, 'max');
    end
    
    %% plot MSE
    %HBO and HbR
    for hh = 1:2
        ttl= [hblab{hh} ' MSE'];
        contour_plots(squeeze(MSE_CCA(:,hh,:,:,:)), ttl,evparams, pOpt, cntno, 'min');
    end
    
    %% plot pvals
    if pvalflag
        %HBO and HbR
        for hh = 1:2
            ttl= [hblab{hh} ' p-values'];
            contour_plots(squeeze(pval_CCA(:,hh,:,:,:)), ttl,evparams, pOpt, cntno, 'min');
        end
    end
    
    %% plot FSCORE
    %HBO and HbR
    for hh = 1:2
        ttl= [hblab{hh} ' FSCORE'];
        contour_plots(squeeze(F_score_CCA(:,hh,:,:,:)), ttl,evparams, pOpt, cntno, 'max');
    end
end


