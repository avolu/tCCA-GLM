% ##### FOLLOWING TWO LINES NEED CHANGE ACCORDING TO USER!
malexflag = 1;
if malexflag
    %Meryem
    path.code = 'C:\Users\mayucel\Documents\PROJECTS\CODES\tCCA-GLM'; addpath(genpath(path.code)); % code directory
    path.dir = 'C:\Users\mayucel\Google Drive\tCCA_GLM_PAPER'; % data directory
    path.save = 'C:\Users\mayucel\Google Drive\tCCA_GLM_PAPER'; % save directory
else
    %Alex
    path.code = 'D:\Office\Research\Software - Scripts\Matlab\Regression tCCA GLM\tCCA-GLM'; addpath(genpath(path.code)); % code directory
    path.dir = 'C:\Users\avolu\Google Drive\tCCA_GLM_PAPER'; % data directory
    path.save = 'C:\Users\avolu\Google Drive\tCCA_GLM_PAPER'; % save directory
end

% #####
filename = 'resting_sim';
set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig))
set(groot,'defaultAxesCreateFcn',@(ax,~)set(ax.Toolbar,'Visible','off'))
sbjfolder = {'Subj33','Subj34','Subj36','Subj37','Subj38','Subj39', 'Subj40', 'Subj41', 'Subj43', 'Subj44','Subj46','Subj47','Subj49','Subj51'};


% Validation parameters
tlags = 0:1:10;
stpsize = 2:2:24;
cthresh = 0.1:0.1:1;

%% load results data from all subjects
CORR_CCA = [];
% #CH x 2(Hbo+HbR) x 2 (cv split) x SBJ x tlag x stepsize x corrthres
for sbj = 1%:numel(sbjfolder)
    res{sbj} = load([path.save '\results_sbj' num2str(sbj) '.mat']);
    
    
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


%% ++++++++++++++++++++++++++++
% THIS IS EXPERIMENTAL AND FOR VALIDATION ATM


%% average across channels
% HERE WE NEED TO REDUCE TO e.g. ONLY TP
CORR_CCA = squeeze(nanmean(CORR_CCA,2));
CORR_SS = squeeze(nanmean(CORR_SS,2));
MSE_CCA = squeeze(nanmean(MSE_CCA,2));
MSE_SS = squeeze(nanmean(MSE_SS,2));
pval_CCA = squeeze(nanmean(pval_CCA,2));
pval_SS = squeeze(nanmean(pval_SS,2));

%% now average across splits
% CORR_CCA = squeeze(nanmean(CORR_CCA,3));
% CORR_SS = squeeze(nanmean(CORR_SS,3));
% MSE_CCA = squeeze(nanmean(MSE_CCA,3));
% MSE_SS = squeeze(nanmean(MSE_SS,3));
% pval_CCA = squeeze(nanmean(pval_CCA,3));
% pval_SS = squeeze(nanmean(pval_SS,3));

%% now average across subjects
% CORR_CCA = squeeze(nanmean(CORR_CCA,1));
% CORR_SS = squeeze(nanmean(CORR_SS,1));
% MSE_CCA = squeeze(nanmean(MSE_CCA,1));
% MSE_SS = squeeze(nanmean(MSE_SS,1));
% pval_CCA = squeeze(nanmean(pval_CCA,1));
% pval_SS = squeeze(nanmean(pval_SS,1));

% for only one subject
CORR_CCA = squeeze(nanmean(CORR_CCA,2));
CORR_SS = squeeze(nanmean(CORR_SS,2));
MSE_CCA = squeeze(nanmean(MSE_CCA,2));
MSE_SS = squeeze(nanmean(MSE_SS,2));
pval_CCA = squeeze(nanmean(pval_CCA,2));
pval_SS = squeeze(nanmean(pval_SS,2));


%% dimensions: HbO/HbR (2) x timelags (11) x stepsize (12) x corr thresh (10)

%% 3D surface plots

x = stpsize;
y = tlags;
z = cthresh;


% plot correlation
figure
[X,Y] = meshgrid(x,y);
surf(X,Y, squeeze(CORR_CCA(1,:,:,1)))
xlabel('stepsize / smpl')
ylabel('time lags / s')
zlabel('HbO correlation')
title('CCA GLM')
figure
[X,Y] = meshgrid(x,y);
surf(X,Y, squeeze(CORR_SS(1,:,:,1)))
xlabel('stepsize / smpl')
ylabel('time lags / s')
zlabel('HbO correlation')
title('SS GLM')

% will be useful, keep for later
%[r,c,v] = ind2sub(size(buf),find(buf == max(buf(:))))



%% visualize # chan
figure
scatter(nump_ss(:), nump_cca(:));
mx=max([nump_ss(:); nump_cca(:)]);
mn=min([nump_ss(:); nump_cca(:)]);
xlim([mn mx])
ylim([mn mx])
hold on
plot([mn mx], [mn mx], 'k')
scatter(nanmean(nump_ss(:)), nanmean(nump_cca(:)), 'xr')
title(['# of significant channels. lag = ' num2str(timelag) 's, ct = ' num2str(param.ct) ' \tau = ' num2str(param.tau)])
xlabel('SS GLM')
ylabel('CCA GLM')

%% visualize pval
figure
scatter(avgp_ss(:), avgp_cca(:));
mx=max([avgp_ss(:); avgp_cca(:)]);
mn=min([avgp_ss(:); avgp_cca(:)]);
xlim([mn mx])
ylim([mn mx])
hold on
plot([mn mx], [mn mx], 'k')
scatter(nanmean(avgp_ss(:)), nanmean(avgp_cca(:)), 'xr')
title(['Average p-val. lag = ' num2str(timelag) 's, ct = ' num2str(param.ct) ' \tau = ' num2str(param.tau)])
xlabel('SS GLM')
ylabel('CCA GLM')


%% plot FP
%% visualize pval
figure
FP_SS = reshape(FP_SS,size(FP_SS,1)*size(FP_SS,2),1);
FP_CCA = reshape(FP_CCA,size(FP_SS,1)*size(FP_SS,2),1);
scatter(FP_SS,FP_CCA);
mx=max([FP_SS(:); FP_CCA(:)]);
mn=min([FP_SS(:); FP_CCA(:)]);
xlim([mn mx])
ylim([mn mx])
hold on
plot([mn mx], [mn mx], 'k')
scatter(nanmean(FP_SS(:)), nanmean(FP_CCA(:)), 'xr')
title(['# of channels with FP, lag = ' num2str(timelag) 's, ct = ' num2str(param.ct) ' \tau = ' num2str(param.tau)])
xlabel('SS GLM')
ylabel('CCA GLM')


%% Briefly eval p val and # det ch improvement
for sbj = 1:numel(sbjfolder)
    for tt = 1:2
        ss_actidx = find(p_SS{sbj,tt});
        cca_actidx = find(p_CCA{sbj,tt});
        
        % number of activated channels
        nump_ss(sbj,tt) = numel(ss_actidx);
        nump_cca(sbj,tt) = numel(cca_actidx);
        % average pval of discovered channels
        format long
        foo = pOxy_SS{sbj,tt};
        if ~isempty(foo)
            avgp_ss(sbj,tt) = nanmean(foo(sbj,ss_actidx));
        end
        foo = pOxy_CCA{sbj,tt};
        if ~isempty(foo)
            avgp_cca(sbj,tt) = nanmean(foo(sbj,cca_actidx));
        end
    end
end

