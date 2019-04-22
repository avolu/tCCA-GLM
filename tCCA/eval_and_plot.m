% ##### FOLLOWING TWO LINES NEED CHANGE ACCORDING TO USER!
<<<<<<< HEAD
malexflag = 1;
=======
malexflag = 0;
>>>>>>> master
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
% Dimensions of output metrics
% #CH x 2(Hbo+HbR) x 2 (cv split) x tlag x stepsize x corrthres
CORR_CCA = [];
<<<<<<< HEAD
% #CH x 2(Hbo+HbR) x 2 (cv split) x SBJ x tlag x stepsize x corrthres
for sbj = 1%:numel(sbjfolder)
=======
for sbj = 1:numel(sbjfolder)
>>>>>>> master
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
CORR_CCA = squeeze(nanmean(CORR_CCA,3));
CORR_SS = squeeze(nanmean(CORR_SS,3));
MSE_CCA = squeeze(nanmean(MSE_CCA,3));
MSE_SS = squeeze(nanmean(MSE_SS,3));
pval_CCA = squeeze(nanmean(pval_CCA,3));
pval_SS = squeeze(nanmean(pval_SS,3));

%% now average across subjects
CORR_CCA = squeeze(nanmean(CORR_CCA,1));
CORR_SS = squeeze(nanmean(CORR_SS,1));
MSE_CCA = squeeze(nanmean(MSE_CCA,1));
MSE_SS = squeeze(nanmean(MSE_SS,1));
pval_CCA = squeeze(nanmean(pval_CCA,1));
pval_SS = squeeze(nanmean(pval_SS,1));

% for only one subject
% CORR_CCA = squeeze(nanmean(CORR_CCA,2));
% CORR_SS = squeeze(nanmean(CORR_SS,2));
% MSE_CCA = squeeze(nanmean(MSE_CCA,2));
% MSE_SS = squeeze(nanmean(MSE_SS,2));
% pval_CCA = squeeze(nanmean(pval_CCA,2));
% pval_SS = squeeze(nanmean(pval_SS,2));


%% dimensions: HbO/HbR (2) x timelags (11) x stepsize (12) x corr thresh (10)

%% 3D surface plots

x = stpsize;
y = tlags;
z = cthresh;


%% plot correlation
figure
[X,Y] = meshgrid(x,y);
surf(X,Y, squeeze(CORR_CCA(1,:,:,1)),'FaceAlpha',0.5)
xlabel('stepsize / smpl')
ylabel('time lags / s')
zlabel('HbO correlation')
title('CCA GLM')
hold on
surf(X,Y, squeeze(CORR_CCA(1,:,:,10)),'FaceAlpha',0.5)
xlabel('stepsize / smpl')
ylabel('time lags / s')
zlabel('HbO correlation')
title('CCA GLM correlation')

[X,Y] = meshgrid(x,y);
figure
for ii=1:9
    subplot(3,3,ii)
    contourf(X,Y, squeeze(CORR_CCA(1,:,:,ii)), 20)
    xlabel('stepsize / smpl')
    ylabel('time lags / s')
    title(['HbO Correlation ctrsh: ' num2str(cthresh(ii))])
    colormap gray
    colorbar
end




%% plot MSE
figure
[X,Y] = meshgrid(x,y);
surf(X,Y, squeeze(MSE_CCA(1,:,:,1)),'FaceAlpha',0.5)
xlabel('stepsize / smpl')
ylabel('time lags / s')
zlabel('HbO MSE')
title('CCA GLM')
hold on
surf(X,Y, squeeze(MSE_CCA(1,:,:,10)),'FaceAlpha',0.5)
xlabel('stepsize / smpl')
ylabel('time lags / s')
zlabel('HbO MSE')
title('CCA GLM MSE')

[X,Y] = meshgrid(x,y);
figure
for ii=1:9
    subplot(3,3,ii)
    contourf(X,Y, squeeze(MSE_CCA(1,:,:,ii)), 20)
    xlabel('stepsize / smpl')
    ylabel('time lags / s')
    title(['HbO MSE ctrsh: ' num2str(cthresh(ii))])
    colormap gray
    colorbar
end


%% plot pvals
figure
[X,Y] = meshgrid(x,y);
surf(X,Y, squeeze(pval_CCA(1,:,:,1)),'FaceAlpha',0.5)
xlabel('stepsize / smpl')
ylabel('time lags / s')
zlabel('HbO pVals')
title('CCA GLM')
hold on
surf(X,Y, squeeze(pval_CCA(1,:,:,10)),'FaceAlpha',0.5)
xlabel('stepsize / smpl')
ylabel('time lags / s')
zlabel('HbO pval')
title('CCA GLM pval')

[X,Y] = meshgrid(x,y);
figure
for ii=1:9
    subplot(3,3,ii)
    contourf(X,Y, squeeze(pval_CCA(1,:,:,ii)), 20)
    xlabel('stepsize / smpl')
    ylabel('time lags / s')
    title(['HbO pVals ctrsh: ' num2str(cthresh(ii))])
    colormap gray
    colorbar
end

% will be useful, keep for later
%[r,c,v] = ind2sub(size(buf),find(buf == max(buf(:))))


