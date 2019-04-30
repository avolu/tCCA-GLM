clear all

%% +++++++++++++++++++++++
%% SCRIPT CONFIGURATION
% +++++++++++++++++++++++
% user: 1 Meryem | 0 Alex
melexflag = 0;
% select which hrf amplitude data: 50 or 100
hrfamp = 50;
% select which MSE type: 1 (average of single trial HRF MSEs), 0: MSE of average HRF 
mseflag = 1;
% Use only true positives for evaluation of metrics
TP_flag = true;
% number of contours in contour plots
cntno = 15;
% plot pvalue results
pvalflag = false;
% plot other metrics
plotmetrics = true;
%% parameters for determining optima
% normalize metrics: 1 X/max | 2 (X-min)/(max-min)
Jparam.nflag = 2;
% smoothing / optimization metrics: 1 mean, 2 median!!! or 3 all channels
Jparam.mtype = 1;
% use mean (1) or median (2) in metric contour plots
mflag = Jparam.mtype;
% Objective function J weights
Jparam.fact.corr = 0;
Jparam.fact.mse =1;
Jparam.fact.pval =0;
Jparam.fact.fscore=1;
Jparam.fact.HbO=1;
Jparam.fact.HbR=1;
% use weighted region of stepsize reg in all directions around evaluation point?
reg.step = 1;%2;
reg.weight =1;%4;
% segmentation approach: threshold for segmentation
Jparam.thresh = 0.7;
% set optimal point per hand to investigate (overwrites opt function
% result), otherwise leave empty
pOptfix =[];
pOptfix = [4 7 6];
%pOptfix = [4 7 8];
%pOptfix = [2 8 9];
plotOptfix = {pOptfix,[3 12 1]};

%% settings to keep in mind
% hrf = 50, Jparam.mtype = 2, fact.corr=1,mse=2,fscore=2 -> Timelag 2, stepsize corr thresh 0.8
% hrf = 50, Jparam.mtype = 3, fact.corr=1,mse=2,fscore=2 -> Timelag 4, stepsize 4, corr thresh 0.5
%                                                           -> vs Timelag 4 stepsize 4 corr thresh 0.2
% pOpt = [5 2 3];


%% Data
% ##### FOLLOWING TWO LINES NEED CHANGE ACCORDING TO USER!
if melexflag
    %Meryem
    path.code = 'C:\Users\mayucel\Documents\PROJECTS\CODES\tCCA-GLM'; addpath(genpath(path.code)); % code directory
    path.dir = 'C:\Users\mayucel\Google Drive\tCCA_GLM_PAPER\FB_RESTING_DATA'; % data directory
    path.save = 'C:\Users\mayucel\Google Drive\tCCA_GLM_PAPER'; % save directory
    path.cvres50 = 'C:\Users\mayucel\Google Drive\tCCA_GLM_PAPER\CV_results_data_50'; % save directory
    path.cvres100 = 'C:\Users\mayucel\Google Drive\tCCA_GLM_PAPER\CV_results_data_100'; % save directory
    path.cvres50stmse = 'C:\Users\mayucel\Google Drive\tCCA_GLM_PAPER\CV_results_data_50_stMSE'; % save directory
    path.cvres100stmse = 'C:\Users\mayucel\Google Drive\tCCA_GLM_PAPER\CV_results_data_100_stMSE'; % save directory
else
    %Alex
    path.code = 'D:\Office\Research\Software - Scripts\Matlab\Regression tCCA GLM\tCCA-GLM'; addpath(genpath(path.code)); % code directory
    path.dir = 'C:\Users\avolu\Google Drive\tCCA_GLM_PAPER\FB_RESTING_DATA'; % data directory
    path.save = 'C:\Users\avolu\Google Drive\tCCA_GLM_PAPER'; % save directory
    path.cvres50 = 'C:\Users\avolu\Google Drive\tCCA_GLM_PAPER\CV_results_data_50'; % save directory
    path.cvres100 = 'C:\Users\avolu\Google Drive\tCCA_GLM_PAPER\CV_results_data_100'; % save directory
    path.cvres50stmse = 'C:\Users\avolu\Google Drive\tCCA_GLM_PAPER\CV_results_data_50_stMSE'; % save directory
    path.cvres100stmse = 'C:\Users\avolu\Google Drive\tCCA_GLM_PAPER\CV_results_data_100_stMSE'; % save directory
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

for hrff=1:2
    hrfamp = hrff*50;
    %% load results data from all subjects
    % Dimensions of output metrics
    % # of sbjs x #CH x 2(Hbo+HbR) x 2 (cv split) x tlag x stepsize x corrthres
    for sbj = 1:numel(sbjfolder)
        switch hrfamp
            case 50
                switch mseflag
                    case 0
                        res{sbj} = load([path.cvres50 '\results_sbj' num2str(sbj) '.mat']);
                    case 1
                        res{sbj} = load([path.cvres50stmse  '\results_sbj' num2str(sbj) '.mat']);
                end
            case 100
                switch mseflag
                    case 0
                        res{sbj} = load([path.cvres100 '\results_sbj' num2str(sbj) '.mat']);
                    case 1
                        res{sbj} = load([path.cvres100stmse  '\results_sbj' num2str(sbj) '.mat']);
                end   
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
    
    %% Calculate True/false positive/negative rates, precision, recall, ...
    tf_errors
    
    
    %% Find Global topology and optimum with objective function, includes segmentation approach
    % calculate objective function output for all input tupel
    fval{hrff} = J_opt(CORR_CCA, MSE_CCA, pval_CCA, F_score_CCA, Jparam ,reg);
    % find optimal parameter set
    [t,s,c] = ind2sub(size(fval{hrff}),find(fval{hrff} == min(fval{hrff}(:))));
    %% overwrite if OPT POINT chosen individually before for further exploration
    if isempty(pOptfix)
        pOpt = [t s c];
        plotOptfix{1} = pOpt;
    else
        pOpt = pOptfix;
        disp('=================================================================')
        disp(['these parameters were chosen manually: ' ...
            num2str(tlags(pOpt(1))) 's, stepsize: ' num2str(stpsize(pOpt(2))) 'smpl, corr threshold: ' num2str(cthresh(pOpt(3)))] )
        disp('=================================================================')
    end
    
    
    disp('=================================================================')
    disp(['these parameters minimize the objective function: timelag: ' ...
        num2str(tlags(t)) 's, stepsize: ' num2str(stpsize(s)) 'smpl, corr threshold: ' num2str(cthresh(c))] )
    disp('=================================================================')
    
    %% Calculate median/mean metrics
    [CORR,MSE,PVAL,FSCORE] = medmean(CORR_CCA, MSE_CCA, pval_CCA, F_score_CCA, mflag);
    
    %% create combined surface plots (depict objective function)
    hblab = {'HbO', 'HbR'};
    
    %% Plot Objective function results
    % normalize fval
    fval{hrff} = (fval{hrff}-min(fval{hrff}(:)))/(max(fval{hrff}(:))-min(fval{hrff}(:)));
    ttl= ['Objective Function, hrf= ' num2str(hrfamp)];
    contour_plots(fval{hrff}, ttl, evparams, pOpt, cntno, 'min');
    
    if plotmetrics
        %% plot correlation
        %HBO and HbR
        for hh = 1:2
            ttl= [hblab{hh} ' Correlation'];
            contour_plots(squeeze(CORR(:,hh,:,:,:)), ttl,evparams, pOpt, cntno, 'max');
        end
        
        %% plot MSE
        %HBO and HbR
        for hh = 1:2
            ttl= [hblab{hh} ' MSE'];
            contour_plots(squeeze(MSE(:,hh,:,:,:)), ttl,evparams, pOpt, cntno, 'min');
        end
        
        %% plot pvals
        if pvalflag
            %HBO and HbR
            for hh = 1:2
                ttl= [hblab{hh} ' p-values'];
                contour_plots(squeeze(PVAL(:,hh,:,:,:)), ttl,evparams, pOpt, cntno, 'min');
            end
        end
        
        %% plot FSCORE
        %HBO and HbR
        for hh = 1:2
            ttl= [hblab{hh} ' FSCORE'];
            contour_plots(squeeze(FSCORE(:,hh,:,:,:)), ttl,evparams, pOpt, cntno, 'max');
        end
    end
    
    %% plot Summary contours for fixed correlation threshold
    ct = pOpt(3);
    [X,Y] = meshgrid(evparams.stpsize,evparams.tlags);
    figure
    dat = {1-fval{hrff}(:,:,ct), squeeze(CORR(:,1,:,:,ct)), squeeze(-MSE(:,1,:,:,ct)), ...
        squeeze(FSCORE(:,1,:,:,ct)), [], ...
        squeeze(CORR(:,2,:,:,ct)), squeeze(-MSE(:,2,:,:,ct)), ...
        squeeze(FSCORE(:,2,:,:,ct))};
    ttl = {'J Opt','CORR HbO','MSE HbO','F HbO', '', 'CORR HbR', 'MSE HbR','F HbR'};
    for dd = 1:numel(dat)
        if ~isempty(dat{dd})
            subplot(2,4,dd)
            climits = [min(dat{dd}(:)) max(dat{dd}(:))];
            contourf(X,Y, dat{dd}, cntno)
            xlabel('stepsize / smpl')
            ylabel('time lags / s')
            title([ttl{dd} ', cthresh: ' num2str(cthresh(ct)), ', hrf=' num2str(hrfamp)])
            colormap hot
            limit = climits(2);
            colorbar
            caxis(climits)
            % mark optimum from objective function
            hold on
            plot(evparams.stpsize(pOpt(1,2)),evparams.tlags(pOpt(1,1)),'diamond','MarkerFaceColor', 'c')
            text(evparams.stpsize(pOpt(1,2)),evparams.tlags(pOpt(1,1)), ['\leftarrow ' num2str(dat{dd}(pOpt(1,1),pOpt(1,2)))])
        end
    end
    
    
    
    %% create scatter plots comparing SS and tCCA
    %append all points
    [CORRcca,MSEcca,PVALcca,FSCOREcca] = medmean(CORR_CCA, MSE_CCA, pval_CCA, F_score_CCA, 3);
    [CORRss,MSEss,PVALss,FSCOREss] = medmean(CORR_SS, MSE_SS, pval_SS, F_score_SS, 3);
    %pOpt = [pOpt(1) pOpt(2) 4]; %(re-)set the optimal parameterset
    figure
    ttl = {'CORR', 'MSE', 'F-SCORE'};
    datss = {squeeze(CORRss(:,:,pOpt(1),pOpt(2),pOpt(3))), squeeze(MSEss(:,:,pOpt(1),pOpt(2),pOpt(3))), squeeze(FSCOREss(:,:,pOpt(1),pOpt(2),pOpt(3)))};
    datcca = {squeeze(CORRcca(:,:,pOpt(1),pOpt(2),pOpt(3))), squeeze(MSEcca(:,:,pOpt(1),pOpt(2),pOpt(3))), squeeze(FSCOREcca(:,:,pOpt(1),pOpt(2),pOpt(3)))};
    ptcol = {'+r', 'xb', '*k'};
    for ff = 1:3
        for hh = 1:2
            axlim = [min([datss{ff}(:,hh); datcca{ff}(:,hh)]) max([datss{ff}(:,hh); datcca{ff}(:,hh)])];
            subplot(2,3,(hh-1)*3+ff)
            hold on
            scatter(squeeze(datss{ff}(:,hh)), squeeze(datcca{ff}(:,hh)), ptcol{hh})
            plot([axlim(1) axlim(2)], [axlim(1) axlim(2)] ,'k')
            scatter(nanmean(squeeze(datss{ff}(:,hh))), nanmean(squeeze(datcca{ff}(:,hh))), ptcol{3})
            % ttest
            [h,p] = ttest(squeeze(datss{ff}(:,hh)),squeeze(datcca{ff}(:,hh)));
            if h
                scatter(nanmean(squeeze(datss{ff}(:,hh))), nanmean(squeeze(datcca{ff}(:,hh))), 'ok')
            end
            title([ttl{ff} ' for Tlag/Ssize/Cthresh: ' num2str(tlags(pOpt(1))) ' / ' num2str(stpsize(pOpt(2))) ' / ' num2str(cthresh(pOpt(3))), ' | p = ' num2str(p)])
            xlim ([axlim(1) axlim(2)])
            ylim ([axlim(1) axlim(2)])
            xlabel('SS GLM')
            ylabel('tCCA GLM')
            grid on
        end
    end
    
    
    %% Plot MSE and F-Score vs corr threshold
    [CORRcca,MSEcca,PVALcca,FSCOREcca] = medmean(CORR_CCA, MSE_CCA, pval_CCA, F_score_CCA, mflag);
    [CORRss,MSEss,PVALss,FSCOREss] = medmean(CORR_SS, MSE_SS, pval_SS, F_score_SS, mflag);
    datss=[];
    datcca=[];
    ptype = {'r','b','--r','--b'};
    pttype = {'or','ob','*r','*b'};
    for pp = 1:2
        datss{pp} = {squeeze(MSEss(:,:,plotOptfix{pp}(1),plotOptfix{pp}(2),:)), squeeze(FSCOREss(:,:,plotOptfix{pp}(1),plotOptfix{pp}(2),:))};
        datcca{pp} = {squeeze(MSEcca(:,:,plotOptfix{pp}(1),plotOptfix{pp}(2),:)), squeeze(FSCOREcca(:,:,plotOptfix{pp}(1),plotOptfix{pp}(2),:))};
    end
    figure
    ylabs={'MSE','F-Score'};
    for mm = 1:2
        subplot(1,2,mm)
        for pp = 1:2
            hold on
            if mm==1
                plot(cthresh, 100*(1-datcca{pp}{mm}(1,:)./datss{pp}{mm}(1,:)), ptype{(2*pp)-1});
                plot(cthresh, 100*(1-datcca{pp}{mm}(2,:)./datss{pp}{mm}(2,:)), ptype{2*pp});
            else
                plot(cthresh, 100*(-1+datcca{pp}{mm}(1,:)./datss{pp}{mm}(1,:)), ptype{(2*pp)-1});
                plot(cthresh, 100*(-1+datcca{pp}{mm}(2,:)./datss{pp}{mm}(2,:)), ptype{2*pp});
            end
            xlabel('Corr threshold')
            ylabel(ylabs{mm})
            title(['Average ' ylabs{mm} ' improvement over SS GLM in % / hrf = ' num2str(hrfamp)])
        end
        plot (cthresh, zeros(numel(cthresh),1), '.k')
        legend(['HbO, @ tlag ' num2str(tlags(plotOptfix{1}(1))) 's, stsize ' num2str(stpsize(plotOptfix{1}(2))*1000/25) 'ms'], ...
            ['HbR, @ tlag ' num2str(tlags(plotOptfix{1}(1))) 's, stsize ' num2str(stpsize(plotOptfix{1}(2))*1000/25) 'ms'], ...
            ['HbO, @ tlag ' num2str(tlags(plotOptfix{2}(1))) 's, stsize ' num2str(stpsize(plotOptfix{2}(2))*1000/25) 'ms'], ...
            ['HbR, @ tlag ' num2str(tlags(plotOptfix{2}(1))) 's, stsize ' num2str(stpsize(plotOptfix{2}(2))*1000/25) 'ms'], ...
            'HbO SS GLM', ...
            'HbR SS GLM', ...
            'Location', 'Best')
    end
    
    
end

plotOptfix = {pOptfix,[4 7 7]};

%% Plot combined Mixed Objective Function contour plot
ttl= 'Sum Obj. Functions hrf=50|100';
% find optimal parameter set
[t,s,c] = ind2sub(size(fval{1}+fval{2}),find(fval{1}+fval{2} == min(fval{1}(:)+fval{2}(:))));
contour_plots((fval{1}+fval{2})/max(fval{1}(:)+fval{2}(:)), ttl, evparams, [t,s,c], cntno, 'min');

disp('=================================================================')
    disp(['these parameters minimize the combined objective functions: timelag: ' ...
        num2str(tlags(t)) 's, stepsize: ' num2str(stpsize(s)) 'smpl, corr threshold: ' num2str(cthresh(c))] )
    disp('=================================================================')


