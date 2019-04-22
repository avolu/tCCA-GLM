clear all

% ##### FOLLOWIG TWO LINES NEED CHANGE ACCRODING TO USER!
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
filename = 'resting_sim';
set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig))
set(groot,'defaultAxesCreateFcn',@(ax,~)set(ax.Toolbar,'Visible','off'))
sbjfolder = {'Subj33','Subj34','Subj36','Subj37','Subj38','Subj39', 'Subj40', 'Subj41', 'Subj43', 'Subj44','Subj46','Subj47','Subj49','Subj51'};


% Validation parameters
tlags = 0:1:10;
stpsize = 2:2:24;
cthresh = 0:0.1:0.9;

%% load results data from all subjects
% Dimensions of output metrics
% #CH x 2(Hbo+HbR) x 2 (cv split) x tlag x stepsize x corrthres
CORR_CCA = [];
for sbj = 1:numel(sbjfolder)
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

% #SBJ x #CH x 2(Hbo+HbR) x 2 (cv split) x tlag x stepsize x corrthres

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


%% dimensions: HbO/HbR (2) x timelags (11) x stepsize (12) x corr thresh (10)

%% 3D surface plots

x = stpsize;
y = tlags;
z = cthresh;
hblab = {'HbO', 'HbR'};

%% plot correlation
%HBO and HbR
for hh = 1:2
    [X,Y] = meshgrid(x,y);
    figure
    climits = [min(min(min(squeeze(CORR_CCA(hh,:,:,:))))) max(max(max(squeeze(CORR_CCA(hh,:,:,:)))))];
    for ii=2:10
        subplot(3,3,ii-1)
        contourf(X,Y, squeeze(CORR_CCA(hh,:,:,ii)), 30)
        xlabel('stepsize / smpl')
        ylabel('time lags / s')
        title([hblab{hh} ' Correlation ctrsh: ' num2str(cthresh(ii))])
        colormap hot
        colorbar
        caxis(climits)
        % mark maxima
        hold on
        buf =  squeeze(CORR_CCA(hh,:,:,ii));
        [r,c] = ind2sub(size(buf),find(buf == max(buf(:))));
        if squeeze(CORR_CCA(hh,r(1),c(1),ii)) == climits(2)
            plot(stpsize(c),tlags(r),'ko','MarkerFaceColor', 'g')
        else
            plot(stpsize(c),tlags(r),'ko','MarkerFaceColor', 'k')
        end
    end
end



%% plot MSE
%HBO & HbR
for hh=1:2
    [X,Y] = meshgrid(x,y);
    figure
    climits = [min(min(min(squeeze(MSE_CCA(hh,:,:,ii))))) max(max(max(squeeze(MSE_CCA(hh,:,:,ii)))))];
    for ii=2:10
        subplot(3,3,ii-1)
        contourf(X,Y, squeeze(MSE_CCA(hh,:,:,ii)), 30)
        xlabel('stepsize / smpl')
        ylabel('time lags / s')
        title([hblab{hh} ' MSE ctrsh: ' num2str(cthresh(ii))])
        colormap(flipud(hot))
        colorbar
        caxis(climits)
        % mark minima
        hold on
        buf =  squeeze(MSE_CCA(hh,:,:,ii));
        [r,c] = ind2sub(size(buf),find(buf == min(buf(:))));
        if squeeze(MSE_CCA(hh,r(1),c(1),ii)) == climits(1)
            plot(stpsize(c),tlags(r),'ko','MarkerFaceColor', 'g')
        else
            plot(stpsize(c),tlags(r),'ko','MarkerFaceColor', 'k')
        end
    end
end

%% plot pvals
%HBO & HbR
for hh=1:2
    [X,Y] = meshgrid(x,y);
    figure
    climits = [min(min(min(squeeze(pval_CCA(hh,:,:,ii))))) max(max(max(squeeze(pval_CCA(hh,:,:,ii)))))];
    for ii=2:10
        subplot(3,3,ii-1)
        contourf(X,Y, squeeze(pval_CCA(hh,:,:,ii)), 20)
        xlabel('stepsize / smpl')
        ylabel('time lags / s')
        title([hblab{hh} ' pVals ctrsh: ' num2str(cthresh(ii))])
        colormap(flipud(hot))
        colorbar
        caxis(climits)
        % mark minima
        hold on
        buf =  squeeze(pval_CCA(hh,:,:,ii));
        [r,c] = ind2sub(size(buf),find(buf == min(buf(:))));
        if squeeze(pval_CCA(hh,r(1),c(1),ii)) == climits(1)
            plot(stpsize(c),tlags(r),'ko','MarkerFaceColor', 'g')
        else
            plot(stpsize(c),tlags(r),'ko','MarkerFaceColor', 'k')
        end
    end
end


% will be useful, keep for later
%[r,c,v] = ind2sub(size(buf),find(buf == max(buf(:))))


%% First try of objective function and finding optimum





%% Plot 3D surface plots
%
% %corr
% figure
% [X,Y] = meshgrid(x,y);
% surf(X,Y, squeeze(CORR_CCA(1,:,:,1)),'FaceAlpha',0.5)
% xlabel('stepsize / smpl')
% ylabel('time lags / s')
% zlabel('HbO correlation')
% title('CCA GLM')
% hold on
% surf(X,Y, squeeze(CORR_CCA(1,:,:,10)),'FaceAlpha',0.5)
% xlabel('stepsize / smpl')
% ylabel('time lags / s')
% zlabel('HbO correlation')
% title('CCA GLM correlation')
% %MSE
% figure
% [X,Y] = meshgrid(x,y);
% surf(X,Y, squeeze(MSE_CCA(1,:,:,1)),'FaceAlpha',0.5)
% xlabel('stepsize / smpl')
% ylabel('time lags / s')
% zlabel('HbO MSE')
% title('CCA GLM')
% hold on
% surf(X,Y, squeeze(MSE_CCA(1,:,:,10)),'FaceAlpha',0.5)
% xlabel('stepsize / smpl')
% ylabel('time lags / s')
% zlabel('HbO MSE')
% title('CCA GLM MSE')
% %pvals
% figure
% [X,Y] = meshgrid(x,y);
% surf(X,Y, squeeze(pval_CCA(1,:,:,1)),'FaceAlpha',0.5)
% xlabel('stepsize / smpl')
% ylabel('time lags / s')
% zlabel('HbO pVals')
% title('CCA GLM')
% hold on
% surf(X,Y, squeeze(pval_CCA(1,:,:,10)),'FaceAlpha',0.5)
% xlabel('stepsize / smpl')
% ylabel('time lags / s')
% zlabel('HbO pval')
% title('CCA GLM pval')