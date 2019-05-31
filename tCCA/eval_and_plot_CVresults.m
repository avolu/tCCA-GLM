clear all

%% +++++++++++++++++++++++
%% SCRIPT CONFIGURATION
% +++++++++++++++++++++++
% user: 1 Meryem | 0 Alex
melexflag = 1;
% select which hrf amplitude data: 1 (20%), 2 (50%) or 3 (100%)
hhh = [1 2 3];
% select which metric type: 1 (average of single trial HRF MSEs), 2: MSE of block average HRF
mmm = [1 2];
% Use only true positives for evaluation of metrics
TP_flag = true;
% number of contours in contour plots
cntno = 15;
% plot pvalue results
pvalflag = false;
% plot other metrics
plotmetrics = false;
% plot corresponding number to local optimum in obj function contour plots
lopttext = false;
% save plots
saveplot = false;
%% parameters for determining optima
% normalize metrics: 1 X/max | 2 (X-min)/(max-min)
Jparam.nflag = 2;
% smoothing / optimization metrics: 1 mean, 2 median!!! or 3 all channels
Jparam.mtype = 1;
% use mean (1) or median (2) in metric contour plots
mflag = Jparam.mtype;
% Objective function J weights
Jparam.fact.corr = 1;
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
pOptfix = [4 8 6];
%pOptfix = [4 7 8];
%pOptfix = [2 8 9];
plotOptfix = {pOptfix,[5 2 1]};

% fixed scatter plot limits (for better visualization, annotate outliers
% per hand afterwards!
sclimflag = true;
sclims = {[-0.2 1], [-0.2 1]; [0 8]*1e-6, [0 4]*1e-6; [0 1], [0 1]};

%% settings to keep in mind
% hrf = 50, Jparam.mtype = 2, fact.corr=1,mse=2,fscore=2 -> Timelag 2, stepsize corr thresh 0.8
% hrf = 50, Jparam.mtype = 3, fact.corr=1,mse=2,fscore=2 -> Timelag 4, stepsize 4, corr thresh 0.5
%                                                           -> vs Timelag 4 stepsize 4 corr thresh 0.2
% pOpt = [5 2 3];

%% get colormaps
cmap_hbo= flipud(othercolor('YlOrRd9'));
cmap_hbr= flipud(othercolor('YlGnBu9'));
cmap_obj= flipud(othercolor('Greys9'));

%% Data
% ##### FOLLOWING TWO LINES NEED CHANGE ACCORDING TO USER!
if melexflag
    %Meryem
    path.code = 'C:\Users\mayucel\Documents\PROJECTS\CODES\tCCA-GLM'; addpath(genpath(path.code)); % code directory
    path.dir = 'C:\Users\mayucel\Google Drive\tCCA_GLM_PAPER\FB_RESTING_DATA'; % data directory
    path.save = 'C:\Users\mayucel\Google Drive\tCCA_GLM_PAPER'; % save directory
    path.cvres20 = 'C:\Users\mayucel\Google Drive\tCCA_GLM_PAPER\CV_results_data_20'; % save directory
    path.cvres50 = 'C:\Users\mayucel\Google Drive\tCCA_GLM_PAPER\CV_results_data_50'; % save directory
    path.cvres100 = 'C:\Users\mayucel\Google Drive\tCCA_GLM_PAPER\CV_results_data_100'; % save directory
    path.cvres20stmse = 'C:\Users\mayucel\Google Drive\tCCA_GLM_PAPER\CV_results_data_20_stMSE'; % save directory
    path.cvres50stmse = 'C:\Users\mayucel\Google Drive\tCCA_GLM_PAPER\CV_results_data_50_stMSE'; % save directory
    path.cvres100stmse = 'C:\Users\mayucel\Google Drive\tCCA_GLM_PAPER\CV_results_data_100_stMSE'; % save directory
    path.savefig = 'C:\Users\mayucel\Google Drive\tCCA_GLM_PAPER\FIGURES\Fig 5-8 contour plots + scatter'
else
    %Alex
    path.code = 'D:\Office\Research\Software - Scripts\Matlab\Regression tCCA GLM\tCCA-GLM'; addpath(genpath(path.code)); % code directory
    path.dir = 'C:\Users\avolu\Google Drive\tCCA_GLM_PAPER\FB_RESTING_DATA'; % data directory
    path.save = 'C:\Users\avolu\Google Drive\tCCA_GLM_PAPER'; % save directory
    path.cvres20 = 'C:\Users\avolu\Google Drive\tCCA_GLM_PAPER\CV_results_data_20'; % save directory
    path.cvres50 = 'C:\Users\avolu\Google Drive\tCCA_GLM_PAPER\CV_results_data_50'; % save directory
    path.cvres100 = 'C:\Users\avolu\Google Drive\tCCA_GLM_PAPER\CV_results_data_100'; % save directory
    path.cvres20stmse = 'C:\Users\avolu\Google Drive\tCCA_GLM_PAPER\CV_results_data_20_stMSE'; % save directory
    path.cvres50stmse = 'C:\Users\avolu\Google Drive\tCCA_GLM_PAPER\CV_results_data_50_stMSE'; % save directory
    path.cvres100stmse = 'C:\Users\avolu\Google Drive\tCCA_GLM_PAPER\CV_results_data_100_stMSE'; % save directory
    path.savefig = 'C:\Users\avolu\Google Drive\tCCA_GLM_PAPER\FIGURES\Fig 5-8 contour plots + scatter'
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

hblab = {'HbO', 'HbR'};
metrttl = {'Single Trial', 'Block Average'};

% xticklabels stepsize in seconds
for tl = 1:numel(evparams.stpsize)
    xtl{tl} = num2str(evparams.stpsize(tl)/25, '%.2g');
    xtl{tl}=strrep(xtl{tl}, '0.', '.');
end

mseffig = figure;

for metr=mmm
    for hrff=hhh
        %% load results data from all subjects
        % Dimensions of output metrics
        % # of sbjs x #CH x 2(Hbo+HbR) x 2 (cv split) x tlag x stepsize x corrthres
        for sbj = 1:numel(sbjfolder)
            switch hrff
                case 1
                    hrfamp = 20;
                    switch metr
                        case 1
                            res{sbj} = load([path.cvres20stmse  '\results_sbj' num2str(sbj) '.mat']);
                        case 2
                            res{sbj} = load([path.cvres20 '\results_sbj' num2str(sbj) '.mat']);
                            
                    end
                case 2
                    hrfamp = 50;
                    switch metr
                        case 1
                            res{sbj} = load([path.cvres50stmse  '\results_sbj' num2str(sbj) '.mat']);
                        case 2
                            res{sbj} = load([path.cvres50 '\results_sbj' num2str(sbj) '.mat']);
                            
                    end
                case 3
                    hrfamp = 100;
                    switch metr
                        case 1
                            res{sbj} = load([path.cvres100stmse  '\results_sbj' num2str(sbj) '.mat']);
                        case 2
                            res{sbj} = load([path.cvres100 '\results_sbj' num2str(sbj) '.mat']);
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
        tf_errors_CV
        
        
        %% Find Global topology and optimum with objective function, includes segmentation approach
        % calculate objective function output for all input tupel
        fval{hrff,metr} = J_opt(CORR_CCA, MSE_CCA, pval_CCA, F_score_CCA, Jparam ,reg);
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
        
        %% Plot Objective function results
        % normalize fval
        fval{hrff,metr} = (fval{hrff,metr}-min(fval{hrff,metr}(:)))/(max(fval{hrff,metr}(:))-min(fval{hrff,metr}(:)));
        ttl= ['Obj. func., hrf= ' num2str(hrfamp) ', ' metrttl{metr}];
        contour_plots(fval{hrff,metr}, ttl, evparams, pOpt, cntno, cmap_obj, 'min', lopttext);
        
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
        dat = {1-fval{hrff,metr}(:,:,ct), squeeze(CORR(:,1,:,:,ct)), squeeze(MSE(:,1,:,:,ct)), ...
            squeeze(FSCORE(:,1,:,:,ct)), [], ...
            squeeze(CORR(:,2,:,:,ct)), squeeze(MSE(:,2,:,:,ct)), ...
            squeeze(FSCORE(:,2,:,:,ct))};
        ttl = {'J Opt','CORR','MSE','F-Score', '', 'CORR', 'MSE', 'F-Score'};
        for dd = 1:numel(dat)
            if ~isempty(dat{dd})
                ax{dd} = subplot(2,4,dd);
                climits = [min(dat{dd}(:)) max(dat{dd}(:))];
                contourf(X,Y, dat{dd}, cntno)
                xlabel('stepsize / s')
                xticks(evparams.stpsize(1:2:end))
                xticklabels(xtl(1:2:end))
                ylabel('time lags / s')
                title([ttl{dd}])
                limit = climits(2);
                if dd ==1
                    colormap(ax{dd},cmap_obj);
                elseif dd<numel(dat)/2+1
                    % MSE?
                    if dd == 3
                        colormap(ax{dd},flipud(cmap_hbo));
                    else
                        colormap(ax{dd},cmap_hbo);
                    end
                else
                    % MSE?
                    if dd == 7
                        colormap(ax{dd},flipud(cmap_hbr));
                    else
                        colormap(ax{dd},cmap_hbr);
                    end
                end
                caxis(climits)
                colorbar
                % mark optimum from objective function
                hold on
                plot(evparams.stpsize(pOpt(1,2)),evparams.tlags(pOpt(1,1)),'diamond','MarkerFaceColor', 'c')
                text(evparams.stpsize(pOpt(1,2)),evparams.tlags(pOpt(1,1)), ['\leftarrow ' num2str(dat{dd}(pOpt(1,1),pOpt(1,2)), '%.2g')]) 
            end
            if dd ==5
                subplot(2,4,dd)
                text(0,0.5,{['HRF: ' num2str(hrfamp) '%']; metrttl{metr}; ['Correlation Threshold: ' num2str(cthresh(ct))]}, 'FontWeight', 'bold'); 
                axis off
            end
        end
        set(gcf, 'Position',  [0,538,1300,458])
        if saveplot
            export_fig([path.savefig '\contours_hrf=' num2str(hrfamp) '%_' metrttl{metr}], '-pdf', '-transparent')
            export_fig([path.savefig '\contours_hrf=' num2str(hrfamp) '%_' metrttl{metr}], '-png', '-transparent', '-r300')
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
                outl = isoutlier(datss{ff}(:,hh));
                axlim = [min([datss{ff}(:,hh); datcca{ff}(:,hh)]) max([datss{ff}(:,hh); datcca{ff}(:,hh)])];
                subplot(2,3,(hh-1)*3+ff)
                hold on
                scatter(squeeze(datss{ff}(:,hh)), squeeze(datcca{ff}(:,hh)), ptcol{hh})
                scatter(nanmean(squeeze(datss{ff}(:,hh))), nanmean(squeeze(datcca{ff}(:,hh))), ptcol{3})
                % ttest
                [h,p] = ttest(squeeze(datss{ff}(:,hh)),squeeze(datcca{ff}(:,hh)));
                if h
                    scatter(nanmean(squeeze(datss{ff}(:,hh))), nanmean(squeeze(datcca{ff}(:,hh))), 'ok')
                end
                title([ttl{ff} ', hrf=' num2str(hrfamp) ', for t/s/c: ' num2str(tlags(pOpt(1))) '/' num2str(stpsize(pOpt(2))) '/' num2str(cthresh(pOpt(3))), ' | p = ' num2str(p, '%0.2g') ', ' metrttl{metr}])
                if sclimflag
                    xlim(sclims{ff,hh})
                    ylim(sclims{ff,hh})
                    plot([sclims{ff,hh}(1) sclims{ff,hh}(2)], [sclims{ff,hh}(1) sclims{ff,hh}(2)] ,'k')
                else
                    xlim ([axlim(1) axlim(2)])
                    ylim ([axlim(1) axlim(2)])
                    plot([axlim(1) axlim(2)], [axlim(1) axlim(2)] ,'k')
                end
                xlabel('SS GLM')
                ylabel('tCCA GLM')
                grid on
            end
        end
        
        
        %% Plot CORR, MSE and F-Score vs corr threshold
        datss = {squeeze(CORRss(:,:,pOpt(1),pOpt(2),:)), squeeze(MSEss(:,:,pOpt(1),pOpt(2),:)), squeeze(FSCOREss(:,:,pOpt(1),pOpt(2),:))};
        datcca = {squeeze(CORRcca(:,:,pOpt(1),pOpt(2),:)), squeeze(MSEcca(:,:,pOpt(1),pOpt(2),:)), squeeze(FSCOREcca(:,:,pOpt(1),pOpt(2),:))};
        ptype = {'r','b','--r','--b'};
        pttype = {'or','ob','*r','*b'};
        figure(mseffig)
        ylabs={'MSE','F-Score'};
        if metr == 1
            midx = 1;
            pofs =0;
        else
            midx = [1 2];
            pofs =1;
        end
        for mm = midx
            subplot(numel(hhh),3,((hrff-1)*numel(hhh))+mm+pofs)
            grid on
            hold on
            switch mm+pofs
                case 1
                    title(['Average ' ylabs{mm} ', hrf = ' num2str(hrfamp) '%, ' metrttl{metr}])
                    ylim([0.5 3]*1e-6);
                case 2
                    title(['Average ' ylabs{mm} ', hrf = ' num2str(hrfamp) '%, ' metrttl{metr}])
                    ylim([0.5 1.5]*1e-7);
                case 3
                    title(['Average ' ylabs{mm} ', hrf = ' num2str(hrfamp) '%'])
                    ylim([0.25 1]);
            end
            errorbar(cthresh, squeeze(nanmean(datcca{mm+1}(:,1,:))), squeeze(nanstd(datcca{mm+1}(:,1,:)))/sqrt(size(datcca{mm+1},1)), 'r');
            errorbar(cthresh, squeeze(nanmean(datcca{mm+1}(:,2,:))), squeeze(nanstd(datcca{mm+1}(:,2,:)))/sqrt(size(datcca{mm+1},1)), 'b');
            plot(cthresh, squeeze(nanmean(datss{mm+1}(:,1,:))), '--r');
            plot(cthresh, squeeze(nanmean(datss{mm+1}(:,2,:))), '--b');
            % mark selected points
            plot(cthresh(pOpt(3)), squeeze(nanmean(datcca{mm+1}(:,1,pOpt(3)))),'ok')
            plot(cthresh(pOpt(3)), squeeze(nanmean(datcca{mm+1}(:,2,pOpt(3)))),'ok')
            xlabel('Corr threshold')
            ylabel(ylabs{mm})
            %plot (ones(10,1)*cthresh(pOpt(3)), zeros(numel(cthresh),1), '.k')
            %             legend(['HbO, @ tlag ' num2str(tlags(plotOptfix{1}(1))) 's, stsize ' num2str(stpsize(plotOptfix{1}(2))*1000/25) 'ms'], ...
            %                 ['HbR, @ tlag ' num2str(tlags(plotOptfix{1}(1))) 's, stsize ' num2str(stpsize(plotOptfix{1}(2))*1000/25) 'ms'], ...
            %                 'HbO, SS GLM', ...
            %                 'HbR, SS GLM', ...
            %                 'Location', 'Best')
        end
        
        
        
        
    end
end

plotOptfix = {pOptfix,[4 7 7]};

%% Plot combined Mixed Objective Function contour plot
ttl= '\Sigma obj. functions';
fvalmxd = zeros(numel(tlags),numel(stpsize),numel(cthresh));
for metr=mmm
    for hrff=hhh
        fvalmxd = fvalmxd + fval{hrff,metr};
    end
end
fvalmxd = fvalmxd/(numel(mmm)+numel(hhh));

% find optimal parameter set
[t,s,c] = ind2sub(size(fvalmxd),find(fvalmxd == min(fvalmxd(:))));
contour_plots((fvalmxd)/max(fvalmxd(:)), ttl, evparams, [t,s,c], cntno, cmap_obj, 'min', lopttext);
set(gcf, 'Position',  [0,538,1300,458])
if saveplot
    export_fig([path.savefig '\sum_optfunct_all.pdf'], '-pdf', '-transparent')
    export_fig([path.savefig '\sum_optfunct_all.pdf'], '-png', '-transparent', '-r300')
end


disp('=================================================================')
disp(['these parameters minimize the combined objective functions: timelag: ' ...
    num2str(tlags(t)) 's, stepsize: ' num2str(stpsize(s)) 'smpl, corr threshold: ' num2str(cthresh(c))] )
disp('=================================================================')


