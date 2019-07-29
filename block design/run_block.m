clear all;

malexflag = 0;
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
filename = 'block_design.nirs';
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
motionflag = false;
% prune
flag_prune = false;
%plot flag
flag_plot = false;


% Optimum CCA parameters (obtained through cross-validation analysis of resting data -see tCCA folder)
timelag = 3;
sts = 16
cthresh = 0.5;

tic;


for sbj = 1:numel(sbjfolder) % loop across subjects
    disp(['subject #' num2str(sbj)]);
    
    %     %% (re-)initialize result matrices
    %     nTrials= NaN(2,numel(tlags),numel(stpsize),numel(cthresh));
    %     DET_SS= NaN(34,2,2,numel(tlags),numel(stpsize),numel(cthresh));
    %     DET_CCA= NaN(34,2,2,numel(tlags),numel(stpsize),numel(cthresh));
    %     pval_SS = NaN(34,2,2,numel(tlags),numel(stpsize),numel(cthresh));
    %     pval_CCA = NaN(34,2,2,numel(tlags),numel(stpsize),numel(cthresh));
    %     MSE_SS = NaN(16,2,2,numel(tlags),numel(stpsize),numel(cthresh));
    %     MSE_CCA = NaN(16,2,2,numel(tlags),numel(stpsize),numel(cthresh));
    %     CORR_SS = NaN(16,2,2,numel(tlags),numel(stpsize),numel(cthresh));
    %     CORR_CCA = NaN(16,2,2,numel(tlags),numel(stpsize),numel(cthresh));
    
    % change to subject directory
    %     cd([path.dir filesep sbjfolder{sbj} filesep]);
    
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
    
    % lpf and get conc data
    dod = hmrIntensity2OD(d);
    dod = hmrBandpassFilt(dod, fq, 0, 0.5);
    dc = hmrOD2Conc(dod, SD, [6 6]);
    
    % resize conc to match the size of d, HbO first half, HbR second half
    foo = [squeeze(dc(:,1,:)),squeeze(dc(:,2,:))];
    
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
    % get d_long and short(and active) now in conc
    d_long = [foo(:,lstLongAct), foo(:,lstLongAct+size(d,2)/2)]; % first half 690 nm; second half 830 nm
    d_short = [foo(:,lstShortAct), foo(:,lstShortAct+size(d,2)/2)]; % first half 690 nm; second half 830 nm
    
    
    %% lowpass filter AUX signals
    AUX = hmrBandpassFilt(AUX, fq, 0, 0.5);
    %% AUX signals
    AUX = [AUX, d_short]; % full AUX = [acc1 acc2 acc3 PPG BP RESP, d_short];
    %% zscore AUX signals
    AUX = zscore(AUX);
    
    
    % remove motion contaminated trials
    if motionflag
        [tIncAuto] = hmrMotionArtifact(dod,fq,SD,ones(size(d,1),1),0.5,1,30,5);
        [s,tRangeStimReject] = enStimRejection(t,s,tIncAuto,ones(size(d,1),1),[-2  10]);
    end
    
    
    %% Perform GLM with SS
    [yavg_ss, yavgstd_ss, tHRF, nTrialsSS, d_ss, yresid_ss, ysum2_ss, beta_ss, yR_ss] = ...
        hmrDeconvHRF_DriftSS(dc, s, t, SD, [], [], [eval_param.HRFmin eval_param.HRFmax], 1, 1, [0.5 0.5], rhoSD_ssThresh, 1, 3, 0);
    
    %% CCA EVAL
    param.tau = sts; %stepwidth for embedding in samples (tune to sample frequency!)
    param.NumOfEmb = ceil(timelag*fq / sts);
    param.ct = cthresh;   % correlation threshold
    %% Perform CCA on training data % AUX = [acc1 acc2 acc3 PPG BP RESP, d_short];
    % use test data of LD channels without synth HRF
    X = d_long;
    flags.shrink = true;
    [REG,  ADD] = rtcca(X,AUX,param,flags);
    
    
    %% Perform GLM with CCA
    [yavg_cca, yavgstd_cca, tHRF, nTrials, d_cca, yresid_cca, ysum2_cca, beta_cca, yR_cca] = ...
        hmrDeconvHRF_DriftSS(dc, s, t, SD, REG, [], [eval_param.HRFmin eval_param.HRFmax], 1, 1, [0.5 0.5], 0, 0, 3, 0);
    
    
    %% list of channels with stimulus
    lst_stim = find(s==1);
    if lst_stim(1) < abs(eval_param.HRFmin) * fq
        lst_stim = lst_stim(2:end);
    end
    if size(s,1) < lst_stim(end) + abs(eval_param.HRFmax) * fq
        lst_stim = lst_stim(1:end-1);
    end
    
    
    %% EVAL / PLOT
    [pval_SS(:,:,sbj), pval_CCA(:,:,sbj)] = ...
        results_eval_block(sbj, d_ss, d_cca, yavg_ss, yavg_cca, tHRF, timelag, sts, cthresh, lst_stim, SD, fq, eval_param, flag_plot, path);
    % Dimensions of output metrics
    %
    
end



%% save data for subject
disp(['saving...'])
save([path.save '\Visual_Paradigm_results' '\p_val_results.mat'], 'pval_SS', 'pval_CCA');

% get total number of sign channels
for i = 1:2 % HbO/R
    for j = 1:size(pval_SS,3) % # of subjects
        lst_sig_SS(i,j) = size(find(pval_SS(:,i,j)<=0.05),1);
        lst_sig_CCA(i,j) = size(find(pval_CCA(:,i,j)<=0.05),1);
    end
end

HbO_Total_Sign_Ch_SS = sum(lst_sig_SS(1,:))

HbO_Total_Sign_Ch_CCA = sum(lst_sig_CCA(1,:))

% HbR_Total_Sign_Ch_SS = sum(lst_sig_SS(2,:))
% 
% HbR_Total_Sign_Ch_CCA = sum(lst_sig_CCA(2,:))


%% Scatter plot HbO/HbR p-value in log scale
Hb = 1; % HbO:1 ; HbR:2
foo1 = squeeze(pval_SS(:,Hb,:));
foo1 = reshape(foo1,1,size(foo1,1)*size(foo1,2));
lst_ss = find(foo1<=0.05);

foo2 = squeeze(pval_CCA(:,Hb,:));
foo2 = reshape(foo2,1,size(foo2,1)*size(foo2,2));
lst_cca = find(foo2<=0.05);

figure;
% common sig channels
scatter(foo1(intersect(lst_ss,lst_cca)),foo2(intersect(lst_ss,lst_cca)));
hold on;

% sign only with SS
scatter(foo1(lst_ss(find(~ismember(lst_ss,intersect(lst_ss,lst_cca))==1))), foo2(lst_ss(find(~ismember(lst_ss,intersect(lst_ss,lst_cca))==1))),'rs');
% sign only with CCA
scatter(foo1(lst_cca(find(~ismember(lst_cca,intersect(lst_ss,lst_cca))==1))), foo2(lst_cca(find(~ismember(lst_cca,intersect(lst_ss,lst_cca))==1))),'gd');
xlabel('pval-SS');
ylabel('pval-CCA');
set(gca,'xscale','log')
set(gca,'yscale','log')
plot([1e-12 1], [1e-12 1],'m');
plot([1e-12 1], [0.05 0.05],'k');
plot([0.05 0.05], [1e-12 1],'k');
% mean value
plot(mean(foo1(union(lst_ss,lst_cca))),mean(foo2(union(lst_ss,lst_cca))),'mx','MarkerSize',10,'LineWidth',2);


legend('\color{blue} Common Sign Channels', '\color{red} Sign. only for GLM with SS','\color{green} Sign. only for CCA');
xlim([1e-6,1])
ylim([1e-6,1])
grid;



%% Scatter plot for number of sign. channels for each subject for SS and CCA methods
figure;
scatter(lst_sig_SS(Hb,:),lst_sig_CCA(Hb,:));
hold;
plot([0 20], [0 20]);
xlabel('# of sign. channels (GLM with SS)');
ylabel('# of sign. channels (GLM with CCA)');
grid

% sign channels
% ss
mean(lst_sig_SS(Hb,:))
std(lst_sig_SS(Hb,:))

% cca
mean(lst_sig_CCA(Hb,:))
std(lst_sig_CCA(Hb,:))

[h,p,c,stats]=ttest(lst_sig_SS(Hb,:),lst_sig_CCA(Hb,:));



% p-values
% ss
mean(foo1(union(lst_ss,lst_cca)))
std(foo1(union(lst_ss,lst_cca)))

% cca
mean(foo2(union(lst_ss,lst_cca)))
std(foo2(union(lst_ss,lst_cca)))

[h,p,c,stats]=ttest(foo1(union(lst_ss,lst_cca)),foo2(union(lst_ss,lst_cca)));
