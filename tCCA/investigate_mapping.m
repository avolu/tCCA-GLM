%% +++++++++++++++++++++++
%% SCRIPT CONFIGURATION
% +++++++++++++++++++++++
% user: 1 Meryem | 0 Alex
melexflag = 0;
% number of contours in contour plots
cntno = 15;
flag_conc = 1; % if 1 CCA inputs are in conc, if 0 CCA inputs are in intensity
flags.pcaf =  [0 0]; % no pca of X or AUX in CCA

%% GLOBAL OPTIMUM TO INVESTIGATE
pOpt = [4 7 6];
% parameters
tlg = pOpt(1);
sts = pOpt(2);
ct = pOpt(3);

tlags = 0:1:10;
stpsize = 2:2:24;
cthresh = 0:0.1:0.9;

%% Data
% ##### FOLLOWING TWO LINES NEED CHANGE ACCORDING TO USER!
if melexflag
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


tic;

disp('=================================================================')
disp(['these parameters were chosen manually: ' ...
    num2str(tlags(pOpt(1))) 's, stepsize: ' num2str(stpsize(pOpt(2))) 'smpl, corr threshold: ' num2str(cthresh(pOpt(3)))] )
disp('=================================================================')

for sbj = 1:numel(sbjfolder) % loop across subjects
    disp(['subject #' num2str(sbj)]);
    
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
        
        param.tau = sts; %stepwidth for embedding in samples (tune to sample frequency!)
        param.NumOfEmb = ceil(tlg*fq / sts);
        
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
        [REG_trn{sbj,tt},  ADD_trn{sbj,tt}] = perf_temp_emb_cca(X,AUX(trnIDX,:),param,flags);
        
        
        %% now use correlation threshold for CCA outside of function to avoid redundant CCA recalculation
        % overwrite: auxiliary cca components that have
        % correlation > ctr
        compindex=find(ADD_trn{sbj,tt}.ccac>ct);
        %overwrite: reduced mapping matrix Av
        ADD_trn{sbj,tt}.Av_red = ADD_trn{sbj,tt}.Av(:,compindex);
        
        
        %% Calculate testig regressors with CCA mapping matrix A from testing
        REG_tst{sbj,tt} = aux_emb*ADD_trn{sbj,tt}.Av_red;
    end
end
toc;



%% re-sort and process mapping matrices




