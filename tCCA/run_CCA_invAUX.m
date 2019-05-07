clear all;

% ##### FOLLOWING TWO LINES NEED CHANGE ACCORDING TO USER!
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
filenames = {'resting_sim_20', 'resting_sim_50', 'resting_sim'} ;
savefold = {'20', '50', '100'; '20_stMSE','50_stMSE','100_stMSE'};
hrffilenames = {'hrf_simdat_20.mat', 'hrf_simdat_50.mat', 'hrf_simdat_100.mat'};
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
flag_plot = false;
% flag for mse/corr for each trial (0 = get mse for average estimated hrf, 1 = get sum of mse for each trial)
flag_trial = 0;
% parameters
tlags = 0:1:10;
stpsize = 2:2:24;
cthresh = 0:0.1:0.9;
% indices of optimal parameterset
pOptfix = [4 7 6];
% hrf type (1/2/3: 20/50/100%)
hrftype =2;

tic;

%% load ground truth hrf
hrf = load([path.code '\sim HRF\' hrffilenames{hrftype}]);

%% create AUX type combination pairs
% indices:  0: None, 1: ACC, 2: PPG, 3: BP, 4: RESP, 5: SS NIRS, 6 ALL
auxinfo.id = [0 1 2 3 4 5];
auxinfo.lab = {'None','ACC','PPG', 'BP', 'RESP', 'SS NIRS', 'ALL'};
auxinfo.auxid = [   0 0 0 0 0 0 1 1 1 1 2 2 2 3 3 4;...
                    1 2 3 4 5 6 2 3 4 5 3 4 5 4 5 5];
aitern = size(auxinfo.auxid,2);

%iteration number
iterno = 1;
totiter = numel(sbjfolder)*2*size(auxinfo.auxid,2);

%% (re-)initialize result matrices
        nTrials= NaN(numel(sbjfolder),34,2,aitern);
        DET_CCA= NaN(numel(sbjfolder),34,2,2,aitern);
        pval_CCA = NaN(numel(sbjfolder),34,2,2,aitern);
        MSE_CCA = NaN(numel(sbjfolder),16,2,2,aitern);
        CORR_CCA = NaN(numel(sbjfolder),16,2,2,aitern);

for sbj = 1:numel(sbjfolder) % loop across subjects
    disp(['subject #' num2str(sbj)]);
    for aa = 1:aitern
        
        
        % change to subject directory
        cd([path.dir filesep sbjfolder{sbj} filesep]);
        
        %% load data
        [fq, t, AUX, d_long, d_short, d0_long, d0_short, d, d0, SD, s, lstLongAct,lstShortAct,lstHrfAdd] = load_nirs(filenames{hrftype},flag_conc);
        
        %% lowpass filter AUX signals
        AUX = hmrBandpassFilt(AUX, fq, 0, 0.5);
        AUXbuf = [AUX, d0_short]; % full AUX
        %% Create AUX signals (combinations)
        % [acc1 acc2 acc3 PPG BP RESP, d_short];
        auxidx = {[], [1:3], [4], [5], [6], [7:6+size(d0_short,2)], [1:size(AUXbuf,2)]};
        AUX = [AUXbuf(:,auxidx{auxinfo.auxid(1,aa)+1}), AUXbuf(:,auxidx{auxinfo.auxid(2,aa)+1})];
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
            
            % Perform GLM with SS
                    [yavg_ss, yavgstd_ss, tHRF, nTrialsSS, d_ss, yresid_ss, ysum2_ss, beta_ss, yR_ss] = ...
                        hmrDeconvHRF_DriftSS(dc{tt}, s(tstIDX,:), t(tstIDX,:), SD, [], [], [eval_param.HRFmin eval_param.HRFmax], 1, 1, [0.5 0.5], rhoSD_ssThresh, 1, 3, 0);
            
            %% CCA EVAL
            %parameters
            timelag = tlags(pOptfix(1));
            sts = stpsize(pOptfix(2));
            ctr = cthresh(pOptfix(3));
            
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
            
            
            disp(['split: ' num2str(tt) ' aux iter: ' num2str(aa)])
            
            %% now use correlation threshold for CCA outside of function to avoid redundant CCA recalculation
            % overwrite: auxiliary cca components that have
            % correlation > ctr
            compindex=find(ADD_trn{tt}.ccac>ctr);
            %overwrite: reduced mapping matrix Av
            ADD_trn{tt}.Av_red = ADD_trn{tt}.Av(:,compindex);
            
            %% Calculate testig regressors with CCA mapping matrix A from testing
            REG_tst = aux_emb*ADD_trn{tt}.Av_red;
            
            %% Perform GLM with CCA
            [yavg_cca, yavgstd_cca, tHRF, nTrials(sbj,tt,aa), d_cca, yresid_cca, ysum2_cca, beta_cca, yR_cca] = ...
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
            [buf1, DET_CCA(sbj,:,:,tt,aa), buf2, ...
                pval_CCA(sbj,:,:,tt,aa), ROCLAB, buf3, MSE_CCA(sbj,:,:,tt,aa), ...
                buf4, CORR_CCA(sbj,:,:,tt,aa)] = ...
                results_eval(sbj, d_ss, d_cca, yavg_ss, yavg_cca, tHRF, timelag, sts, ctr, lst_stim, SD, fq, lstHrfAdd, eval_param, flag_plot, path, hrf, flag_trial, nTrials(sbj,tt,aa));
            % Dimensions of output metrics
            % # Subjects x #CH x 2(Hbo+HbR) x 2 (cv split) x AUX combinations
            
            % display iterno
            disp(['iter #' num2str(iterno) ', sbj ' num2str(sbj) ', ' num2str(ceil(1000*iterno/(totiter))/10) '% done'])
            iterno = iterno+1;
            
        end
    end
    % clear vars
    clear vars AUX d d0 d_long d0_long d_short d0_short t s REG_trn ADD_trn
end

%% save data
save([path.save '\CV_AUX_contributions_' savefold{flag_trial+1,hrftype} '\results_all_sbj.mat'], 'DET_CCA', 'pval_CCA', 'ROCLAB', 'MSE_CCA', 'CORR_CCA', 'nTrials', 'pOptfix', 'auxinfo');

toc;






