clear all

% ##### FOLLOWING TWO LINES NEED CHANGE ACCORDING TO USER!
malexflag = 0;
if malexflag
    %Meryem
    path.code = 'C:\Users\mayucel\Documents\PROJECTS\CODES\tCCA-GLM'; addpath(genpath(path.code)); % code directory
    path.dir = 'C:\Users\mayucel\Google Drive\tCCA_GLM_PAPER\FB_RESTING_DATA'; % data directory
    path.save = 'C:\Users\mayucel\Google Drive\tCCA_GLM_PAPER'; % save directory
    path.cvres = 'C:\Users\mayucel\Google Drive\tCCA_GLM_PAPER\CV results data'; % save directory
else
    %Alex
    path.code = 'D:\Office\Research\Software - Scripts\Matlab\Regression tCCA GLM\tCCA-GLM'; addpath(genpath(path.code)); % code directory
    path.dir = 'C:\Users\avolu\Google Drive\tCCA_GLM_PAPER\FB_RESTING_DATA'; % data directory
    path.save = 'C:\Users\avolu\Google Drive\tCCA_GLM_PAPER'; % save directory
    path.cvres = 'C:\Users\avolu\Google Drive\tCCA_GLM_PAPER\CV results data'; % save directory
end

% #####
filename = 'resting_sim';
set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig))
set(groot,'defaultAxesCreateFcn',@(ax,~)set(ax.Toolbar,'Visible','off'))
sbjfolder = {'Subj33','Subj34','Subj36','Subj37','Subj38','Subj39', 'Subj40', 'Subj41', 'Subj43', 'Subj44','Subj46','Subj47','Subj49','Subj51'};

% flags
TP_flag = true;

% Validation parameters
tlags = 0:1:10;
stpsize = 2:2:24;
cthresh = 0:0.1:0.9;

%% load results data from all subjects
% Dimensions of output metrics
% # of sbjs x #CH x 2(Hbo+HbR) x 2 (cv split) x tlag x stepsize x corrthres

CORR_CCA = [];
for sbj = 1:numel(sbjfolder)
    res{sbj} = load([path.cvres '\results_sbj' num2str(sbj) '.mat']);
    
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

% True Positive *Rate* and False Positive *Rate*
TPR_SS = Ch_TP_SS./(Ch_TP_SS + Ch_FN_SS);
FPR_SS = Ch_FP_SS./(Ch_FP_SS + Ch_TN_SS);

fPR_CCA = Ch_TP_CCA./(Ch_TP_CCA + Ch_FN_CCA);
FPR_CCA = Ch_FP_CCA./(Ch_FP_CCA + Ch_TN_CCA);


% get F-score
% SS & CCA
Precision_SS = Ch_TP_SS ./(Ch_TP_SS + Ch_FP_SS);
Recall_SS = Ch_TP_SS ./(Ch_TP_SS + Ch_FN_SS);
F_score_SS = 2 * (Precision_SS .* Recall_SS)./(Precision_SS + Recall_SS);

Precision_CCA = Ch_TP_CCA ./(Ch_TP_CCA + Ch_FP_CCA);
Recall_CCA = Ch_TP_CCA ./(Ch_TP_CCA + Ch_FN_CCA);
F_score_CCA = 2 * (Precision_CCA .* Recall_CCA)./(Precision_CCA + Recall_CCA);

% reshape all to "# of sbjs x 2(Hbo+HbR) x 2 (cv split) x tlag x stepsize x corrthres
Ch_TP_SS = reshape(Ch_TP_SS, size(DET_SS,1), size(DET_SS,3), size(DET_SS,4), size(DET_SS,5), size(DET_SS,6), size(DET_SS,7));
Ch_FP_SS = reshape(Ch_FP_SS, size(DET_SS,1), size(DET_SS,3), size(DET_SS,4), size(DET_SS,5), size(DET_SS,6), size(DET_SS,7));
Ch_FN_SS = reshape(Ch_FN_SS, size(DET_SS,1), size(DET_SS,3), size(DET_SS,4), size(DET_SS,5), size(DET_SS,6), size(DET_SS,7));
Ch_TN_SS = reshape(Ch_TN_SS, size(DET_SS,1), size(DET_SS,3), size(DET_SS,4), size(DET_SS,5), size(DET_SS,6), size(DET_SS,7));
F_score_SS = reshape(F_score_SS, size(DET_SS,1), size(DET_SS,3), size(DET_SS,4), size(DET_SS,5), size(DET_SS,6), size(DET_SS,7));
TPR_SS = reshape(TPR_SS, size(DET_SS,1), size(DET_SS,3), size(DET_SS,4), size(DET_SS,5), size(DET_SS,6), size(DET_SS,7));
FPR_SS = reshape(FPR_SS, size(DET_SS,1), size(DET_SS,3), size(DET_SS,4), size(DET_SS,5), size(DET_SS,6), size(DET_SS,7));

Ch_TP_CCA = reshape(Ch_TP_CCA, size(DET_CCA,1), size(DET_CCA,3), size(DET_CCA,4), size(DET_CCA,5), size(DET_CCA,6), size(DET_CCA,7));
Ch_FP_CCA = reshape(Ch_FP_CCA, size(DET_CCA,1), size(DET_CCA,3), size(DET_CCA,4), size(DET_CCA,5), size(DET_CCA,6), size(DET_CCA,7));
Ch_FN_CCA = reshape(Ch_FN_CCA, size(DET_CCA,1), size(DET_CCA,3), size(DET_CCA,4), size(DET_CCA,5), size(DET_CCA,6), size(DET_CCA,7));
Ch_TN_CCA = reshape(Ch_TN_CCA, size(DET_CCA,1), size(DET_CCA,3), size(DET_CCA,4), size(DET_CCA,5), size(DET_CCA,6), size(DET_CCA,7));
F_score_CCA = reshape(F_score_CCA, size(DET_CCA,1), size(DET_CCA,3), size(DET_CCA,4), size(DET_CCA,5), size(DET_CCA,6), size(DET_CCA,7));
fPR_CCA = reshape(fPR_CCA, size(DET_SS,1), size(DET_SS,3), size(DET_SS,4), size(DET_SS,5), size(DET_SS,6), size(DET_SS,7));
FPR_CCA = reshape(FPR_CCA, size(DET_SS,1), size(DET_SS,3), size(DET_SS,4), size(DET_SS,5), size(DET_SS,6), size(DET_SS,7));



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
F_score_CCA = squeeze(nanmean(F_score_CCA,3));
F_score_SS = squeeze(nanmean(F_score_SS,3));
Ch_FP_CCA = squeeze(nanmean(Ch_FP_CCA,3));



%% now average across subjects
CORR_CCA = squeeze(nanmean(CORR_CCA,1));
CORR_SS = squeeze(nanmean(CORR_SS,1));
MSE_CCA = squeeze(nanmean(MSE_CCA,1));
MSE_SS = squeeze(nanmean(MSE_SS,1));
pval_CCA = squeeze(nanmean(pval_CCA,1));
pval_SS = squeeze(nanmean(pval_SS,1));
F_score_CCA = squeeze(nanmean(F_score_CCA,1));
F_score_SS = squeeze(nanmean(F_score_SS,1));
Ch_FP_CCA = squeeze(nanmean(Ch_FP_CCA,1));


%% dimensions: HbO/HbR (2) x timelags (11) x stepsize (12) x corr thresh (10)
x = stpsize;
y = tlags;
z = cthresh;
hblab = {'HbO', 'HbR'};

%% Global optimum using clustering approach: find parameter set(s) that are in each optimal cluster
% and find overlaps
% normalize inputs
for hh=1:2
    CORR(hh,:,:,:) = CORR_CCA(hh,:,:,:)./max(CORR_CCA(:));
    MSE(hh,:,:,:) = MSE_CCA(hh,:,:,:)./max(MSE_CCA(:));
    PVAL(hh,:,:,:) = pval_CCA(hh,:,:,:)./max(pval_CCA(:));
    FSCORE(hh,:,:,:) = F_score_CCA(hh,:,:,:)./max(F_score_CCA(:));
end
% HbO & HbR
for hh=1:2
    k = 30; %start number of clusters
    R{hh} = [];
    while isempty(R{hh}) && k >0
        [L,centers] = imsegkmeans3(single(squeeze(CORR(hh,:,:,:))),k);
        [m,i] = max(centers);
        [t,s,c] = ind2sub(size(L),find(L==i));
        C = [t s c];
        
        [L,centers] = imsegkmeans3(single(squeeze(MSE(hh,:,:,:))),k);
        [m,i] = min(centers);
        [t,s,c] = ind2sub(size(L),find(L==i));
        M = [t s c];
        
        [L,centers] = imsegkmeans3(single(squeeze(PVAL(hh,:,:,:))),k);
        [m,i] = min(centers);
        [t,s,c] = ind2sub(size(L),find(L==i));
        P = [t s c];
        
        [L,centers] = imsegkmeans3(single(squeeze(FSCORE(hh,:,:,:))),k);
        [m,i] = max(centers);
        [t,s,c] = ind2sub(size(L),find(L==i));
        F = [t s c];
        
        R1 = intersect(C,M,'rows');
        R2 = intersect(P,F,'rows');
        R{hh} = intersect(R1,R2,'rows');
        if isempty(R{hh})
            k = k-1;
        else
            clust(hh) = k;
        end
    end
end
buf = intersect(R{1}, R{2}, 'rows');
R=[];
R=buf;

disp('=================================================================')
disp(['Results of cluster analysis, Optimal parameters for HbO and HbR (' num2str(clust(1)) '/' num2str(clust(2)) ' clusters)'])
for ii = 1:size(R,1)
    disp(['Timelag: ' num2str(tlags(R(ii,1))) 's, Stepsize: ' num2str(stpsize(R(ii,2))) 'smpls, Corr Thresh: ' num2str(cthresh(R(ii,3)))])
    disp(['Corresponding CORR/MSE/PVAL/FSCORE: ' num2str(CORR_CCA(hh,R(ii,1),R(ii,2),R(ii,3))) '/ ' ...
        num2str(MSE_CCA(hh,R(ii,1),R(ii,2),R(ii,3))) '/ ' ...
        num2str(pval_CCA(hh,R(ii,1),R(ii,2),R(ii,3))) '/ ' ...
        num2str(F_score_CCA(hh,R(ii,1),R(ii,2),R(ii,3)))])
end
disp('=================================================================')


%% Approach to find a global optimum: use objective function
% normalize inputs
for hh=1:2
    CORR(hh,:,:,:) = CORR_CCA(hh,:,:,:)./max(CORR_CCA(:));
    MSE(hh,:,:,:) = MSE_CCA(hh,:,:,:)./max(MSE_CCA(:));
    PVAL(hh,:,:,:) = pval_CCA(hh,:,:,:)./max(pval_CCA(:));
    FSCORE(hh,:,:,:) = F_score_CCA(hh,:,:,:)./max(F_score_CCA(:));
end
% fact: struct with factors (weights) for adapting the objective function J
fact.corr = 1;
fact.mse =2;
fact.pval =1;
fact.fscore=2;
fact.HbO=1;
fact.HbR=1;
% calculate objective function output for all input tupel
for tt = 1:11
    for ss = 1:12
        for cc = 1:10
            xx=[tt ss cc];
            fval(tt,ss,cc) = J_opt(xx, CORR, MSE, PVAL, FSCORE, fact);
        end
    end
end
% find optimal parameter set
[t,s,c] = ind2sub(size(fval),find(fval == min(fval(:))));
pOpt = [t s c];
disp('=================================================================')
disp(['these parameters minimize the objective function: timelag: ' ...
    num2str(tlags(t)) 's, stepsize: ' num2str(stpsize(s)) 'smpl, corr threshold: ' num2str(cthresh(c))] )
disp('=================================================================')


%% create combined surface plots (depict objective function)
[X,Y] = meshgrid(x,y);
figure
climits = [min(fval(:)) max(fval(:))];
for ii=1:10
    subplot(3,4,ii)
    contourf(X,Y, squeeze(fval(:,:,ii)), 20)
    xlabel('stepsize / smpl')
    ylabel('time lags / s')
    title(['Combined (J), ctrsh: ' num2str(cthresh(ii))])
    colormap(flipud(hot))
    colorbar
    caxis(climits)
    % mark local optima
    hold on
    buf =  squeeze(fval(:,:,ii));
    [r,c] = ind2sub(size(buf),find(buf == min(buf(:))));
    if squeeze(fval(r(1),c(1),ii)) == climits(1)
        plot(stpsize(c),tlags(r),'ko','MarkerFaceColor', 'g')
    else
        plot(stpsize(c),tlags(r),'ko','MarkerFaceColor', 'k')
    end
    % mark optima from cluster analysis
    ridx = find(R(:,3)==ii);
    if ~isempty(ridx)
        for rr = 1:numel(ridx)
            plot(stpsize(R(ridx(rr),2)),tlags(R(ridx(rr),1)),'square','MarkerFaceColor', 'b')
        end
    end
    % mark optimum from objective function
    if ii == pOpt(1,3)
        plot(stpsize(pOpt(1,2)),tlags(pOpt(1,1)),'diamond','MarkerFaceColor', 'c')
    end
end


%% plot correlation
%HBO and HbR
for hh = 1:2
    [X,Y] = meshgrid(x,y);
    figure
    climits = [min(min(min(squeeze(CORR_CCA(hh,:,:,:))))) max(max(max(squeeze(CORR_CCA(hh,:,:,:)))))];
    for ii=1:10
        subplot(3,4,ii)
        contourf(X,Y, squeeze(CORR_CCA(hh,:,:,ii)), 30)
        xlabel('stepsize / smpl')
        ylabel('time lags / s')
        title([hblab{hh} ' Correlation ctrsh: ' num2str(cthresh(ii))])
        colormap hot
        colorbar
        caxis(climits)
        % mark local optima
        hold on
        buf =  squeeze(CORR_CCA(hh,:,:,ii));
        [r,c] = ind2sub(size(buf),find(buf == max(buf(:))));
        if squeeze(CORR_CCA(hh,r(1),c(1),ii)) == climits(2)
            plot(stpsize(c),tlags(r),'ko','MarkerFaceColor', 'g')
        else
            plot(stpsize(c),tlags(r),'ko','MarkerFaceColor', 'k')
        end
        % mark optima from cluster analysis
        ridx = find(R(:,3)==ii);
        if ~isempty(ridx)
            for rr = 1:numel(ridx)
                plot(stpsize(R(ridx(rr),2)),tlags(R(ridx(rr),1)),'square','MarkerFaceColor', 'b')
            end
        end
        % mark optimum from objective function
        if ii == pOpt(1,3)
            plot(stpsize(pOpt(1,2)),tlags(pOpt(1,1)),'diamond','MarkerFaceColor', 'c')
        end
    end
end



%% plot MSE
%HBO & HbR
for hh=1:2
    [X,Y] = meshgrid(x,y);
    figure
    climits = [min(min(min(squeeze(MSE_CCA(hh,:,:,ii))))) max(max(max(squeeze(MSE_CCA(hh,:,:,ii)))))];
    for ii=1:10
        subplot(3,4,ii)
        contourf(X,Y, squeeze(MSE_CCA(hh,:,:,ii)), 30)
        xlabel('stepsize / smpl')
        ylabel('time lags / s')
        title([hblab{hh} ' MSE ctrsh: ' num2str(cthresh(ii))])
        colormap(flipud(hot))
        colorbar
        caxis(climits)
        % mark local optima
        hold on
        buf =  squeeze(MSE_CCA(hh,:,:,ii));
        [r,c] = ind2sub(size(buf),find(buf == min(buf(:))));
        if squeeze(MSE_CCA(hh,r(1),c(1),ii)) == climits(1)
            plot(stpsize(c),tlags(r),'ko','MarkerFaceColor', 'g')
        else
            plot(stpsize(c),tlags(r),'ko','MarkerFaceColor', 'k')
        end
        % mark optima from cluster analysis
        ridx = find(R(:,3)==ii);
        if ~isempty(ridx)
            for rr = 1:numel(ridx)
                plot(stpsize(R(ridx(rr),2)),tlags(R(ridx(rr),1)),'square','MarkerFaceColor', 'b')
            end
        end
        % mark optimum from objective function
        if ii == pOpt(1,3)
            plot(stpsize(pOpt(1,2)),tlags(pOpt(1,1)),'diamond','MarkerFaceColor', 'c')
        end
    end
end

%% plot pvals
%HBO & HbR
for hh=1:2
    [X,Y] = meshgrid(x,y);
    figure
    climits = [min(min(min(squeeze(pval_CCA(hh,:,:,ii))))) max(max(max(squeeze(pval_CCA(hh,:,:,ii)))))];
    for ii=1:10
        subplot(3,4,ii)
        contourf(X,Y, squeeze(pval_CCA(hh,:,:,ii)), 20)
        xlabel('stepsize / smpl')
        ylabel('time lags / s')
        title([hblab{hh} ' pVals ctrsh: ' num2str(cthresh(ii))])
        colormap(flipud(hot))
        colorbar
        caxis(climits)
        % mark local optima
        hold on
        buf =  squeeze(pval_CCA(hh,:,:,ii));
        [r,c] = ind2sub(size(buf),find(buf == min(buf(:))));
        if squeeze(pval_CCA(hh,r(1),c(1),ii)) == climits(1)
            plot(stpsize(c),tlags(r),'ko','MarkerFaceColor', 'g')
        else
            plot(stpsize(c),tlags(r),'ko','MarkerFaceColor', 'k')
        end
        % mark optima from cluster analysis
        ridx = find(R(:,3)==ii);
        if ~isempty(ridx)
            for rr = 1:numel(ridx)
                plot(stpsize(R(ridx(rr),2)),tlags(R(ridx(rr),1)),'square','MarkerFaceColor', 'b')
            end
        end
        % mark optimum from objective function
        if ii == pOpt(1,3)
            plot(stpsize(pOpt(1,2)),tlags(pOpt(1,1)),'diamond','MarkerFaceColor', 'c')
        end
    end
end


%% plot FSCORE
%HBO and HbR
for hh = 1:2
    [X,Y] = meshgrid(x,y);
    figure
    climits = [min(min(min(squeeze(F_score_CCA(hh,:,:,:))))) max(max(max(squeeze(F_score_CCA(hh,:,:,:)))))];
    for ii=1:10
        subplot(3,4,ii)
        contourf(X,Y, squeeze(F_score_CCA(hh,:,:,ii)), 30)
        xlabel('stepsize / smpl')
        ylabel('time lags / s')
        title([hblab{hh} ' FSCORE ctrsh: ' num2str(cthresh(ii))])
        colormap hot
        colorbar
        caxis(climits)
        % mark local optima
        hold on
        buf =  squeeze(F_score_CCA(hh,:,:,ii));
        [r,c] = ind2sub(size(buf),find(buf == max(buf(:))));
        if squeeze(F_score_CCA(hh,r(1),c(1),ii)) == climits(2)
            plot(stpsize(c),tlags(r),'ko','MarkerFaceColor', 'g')
        else
            plot(stpsize(c),tlags(r),'ko','MarkerFaceColor', 'k')
        end
        % mark optima from cluster analysis
        ridx = find(R(:,3)==ii);
        if ~isempty(ridx)
            for rr = 1:numel(ridx)
                plot(stpsize(R(ridx(rr),2)),tlags(R(ridx(rr),1)),'square','MarkerFaceColor', 'b')
            end
        end
        % mark optimum from objective function
        if ii == pOpt(1,3)
            plot(stpsize(pOpt(1,2)),tlags(pOpt(1,1)),'diamond','MarkerFaceColor', 'c')
        end
    end
end



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