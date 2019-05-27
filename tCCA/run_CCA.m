clear all;
 
% ##### FOLLOWING TWO LINES NEED CHANGE ACCORDING TO USER!
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
filename = 'resting_sim_50';
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
 
% Validation parameters
tlags = 3;%0:1:10;
stpsize = 16;%2:2:24;
cthresh = 0.5;%0:0.1:0.9;
 
tlidx =0;
stpidx =0;
ctidx =0;
 
tic;
 
%% load ground truth hrf
hrf = load([path.code '\sim HRF\hrf_simdat_50.mat']);
 
%iteration number
iterno = 1;
totiter = numel(sbjfolder)*2*numel(tlags)*numel(stpsize)*numel(cthresh);
 
for sbj = 4%1:numel(sbjfolder) % loop across subjects
    disp(['subject #' num2str(sbj)]);
    
    %% (re-)initialize result matrices
    nTrials= NaN(2,numel(tlags),numel(stpsize),numel(cthresh));
    DET_SS= NaN(34,2,2,numel(tlags),numel(stpsize),numel(cthresh));
    DET_CCA= NaN(34,2,2,numel(tlags),numel(stpsize),numel(cthresh));
    pval_SS = NaN(34,2,2,numel(tlags),numel(stpsize),numel(cthresh));
    pval_CCA = NaN(34,2,2,numel(tlags),numel(stpsize),numel(cthresh));
    MSE_SS = NaN(16,2,2,numel(tlags),numel(stpsize),numel(cthresh));
    MSE_CCA = NaN(16,2,2,numel(tlags),numel(stpsize),numel(cthresh));
    CORR_SS = NaN(16,2,2,numel(tlags),numel(stpsize),numel(cthresh));
    CORR_CCA = NaN(16,2,2,numel(tlags),numel(stpsize),numel(cthresh));
    
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
        for tl = tlags %loop across timelags
            timelag = tl;
            tlidx = tlidx+1;
            
            for sts = stpsize  %loop across stepsizes
                stpidx = stpidx+1;
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
                
                for ctr = cthresh %loop across correlation thresholds
                    ctidx = ctidx+1;
                    disp(['split: ' num2str(tt) ', tlag: ' num2str(tl) ', stsize: ' num2str(sts) ', ctrhesh: ' num2str(ctr)])
                    
                    %% now use correlation threshold for CCA outside of function to avoid redundant CCA recalculation
                    % overwrite: auxiliary cca components that have
                    % correlation > ctr
                    compindex=find(ADD_trn{tt}.ccac>ctr);
                    %overwrite: reduced mapping matrix Av
                    ADD_trn{tt}.Av_red = ADD_trn{tt}.Av(:,compindex);
                    
                    %% Calculate testig regressors with CCA mapping matrix A from testing
                    REG_tst = aux_emb*ADD_trn{tt}.Av_red;
                    
                    %% Perform GLM with CCA
                    [yavg_cca, yavgstd_cca, tHRF, nTrials(tt,tlidx,stpidx,ctidx), d_cca, yresid_cca, ysum2_cca, beta_cca, yR_cca] = ...
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
                    [DET_SS(:,:,tt,tlidx,stpidx,ctidx), DET_CCA(:,:,tt,tlidx,stpidx,ctidx), pval_SS(:,:,tt,tlidx,stpidx,ctidx), ...
                        pval_CCA(:,:,tt,tlidx,stpidx,ctidx), ROCLAB, MSE_SS(:,:,tt,tlidx,stpidx,ctidx), MSE_CCA(:,:,tt,tlidx,stpidx,ctidx), ...
                        CORR_SS(:,:,tt,tlidx,stpidx,ctidx), CORR_CCA(:,:,tt,tlidx,stpidx,ctidx)] = ...
                        results_eval(sbj, d_ss, d_cca, yavg_ss, yavg_cca, tHRF, timelag, sts, ctr, lst_stim, SD, fq, lstHrfAdd, lstLongAct, eval_param, flag_plot, path, hrf, flag_trial, nTrials(tt,tlidx,stpidx,ctidx));
                    % Dimensions of output metrics
                    % #CH x 2(Hbo+HbR) x 2 (cv split) x tlag x stepsize x corrthres
                    % old:  #CH x 2(Hbo+HbR) x 2 (cv split) x SBJ x tlag x stepsize x corrthres
                    
                    % display iterno
                    disp(['iter #' num2str(iterno) ', sbj ' num2str(sbj) ', ' num2str(ceil(1000*iterno/(totiter))/10) '% done'])
                    iterno = iterno+1;
                end
                % reset counter
                ctidx =0;
            end
            % reset counter
            stpidx =0;
        end
        % reset counter
        tlidx =0;
    end
    %% save data for subject
    disp(['saving sbj ' num2str(sbj) '...'])
    save([path.save '\CV_results_data_20' '\results_sbj' num2str(sbj) '.mat'], 'DET_SS', 'DET_CCA', 'pval_SS', 'pval_CCA', 'ROCLAB', 'MSE_SS', 'MSE_CCA', 'CORR_SS', 'CORR_CCA', 'nTrials');
    % clear vars
    clear vars AUX d d0 d_long d0_long d_short d0_short t s REG_trn ADD_trn
    
end
 
 
toc;
 
 
 
 
 
 

