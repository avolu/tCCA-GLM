%clear all;

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
%% simulated data file names
filename = 'resting_sim_50';
%% load ground truth hrf
hrf = load([path.code '\sim HRF\hrf_simdat_50.mat']);

set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig))
set(groot,'defaultAxesCreateFcn',@(ax,~)set(ax.Toolbar,'Visible','off'))
sbjfolder = {'Subj33','Subj34','Subj36','Subj37','Subj38','Subj39', 'Subj40', 'Subj41', 'Subj43', 'Subj44','Subj46','Subj47','Subj49','Subj51'};


%% Options/Parameter Settings
% CCA parameters
flags.pcaf =  [0 0]; % no pca of X or AUX
flag_conc = 1; % if 1 CCA inputs are in conc, if 0 CCA inputs are in intensity
% Validation parameters
tlags = 3;%0:1:10;
stpsize = 16;%2:2:24;
cthresh = 0.5;%0:0.1:0.9;

for sbj = 1%:numel(sbjfolder) % loop across subjects
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
    for tt = 1%:2
        tstIDX = spltIDX{trntst{tt}(1)};
        trnIDX = spltIDX{trntst{tt}(2)};
        
        param.tau = stpsize; %stepwidth for embedding in samples (tune to sample frequency!)
        param.NumOfEmb = ceil(tlags*fq / stpsize);
        
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
        
        %% %% NEW!!! ZSCORE THE DATA
        X = zscore(d0_long(trnIDX,:));
        % use the test data of aux
        Y = zscore(aux_emb);
        
        %% This is where the following function would run
        [REG_trn{tt},  ADD_trn{tt}] = perf_temp_emb_cca(X,AUX(trnIDX,:),param,flags);
        %% This is the new function with shrinkage
        %% we formulate the generalised eigenvalue equation to investigate the generalized eigenvalue probolem and eigenspectrum
        flags.shrink = false;
        [REG_trn2{tt},  ADD_trn2{tt}] = rtcca(X,AUX(trnIDX,:),param,flags);
        
        %% plot example CCA component of both modalities for both methods
        % black: matlab cca results
        % green: regularized own cca results
        comp=1;
        % find common sign
        s = sign(corr(ADD_trn2{tt}.U(:,comp),ADD_trn{tt}.U(:,comp)));
        figure
        plot(ADD_trn2{tt}.U(:,comp),'g')
        hold on
        plot(ADD_trn2{tt}.V(:,comp),'--g')
        plot(s*ADD_trn{tt}.U(:,comp),'k')
        plot(s*ADD_trn{tt}.V(:,comp),'--k')
        title(["CCA component " num2str(comp) " of both modalities. Regularized: Green."])
        
        % plot canonical correlation coefficients
        figure
        plot(ADD_trn{tt}.ccac,'k')
        hold on
        plot(comp,ADD_trn{tt}.ccac(comp),'ok')
        plot(ADD_trn2{tt}.ccac,'g')
        plot(comp,ADD_trn2{tt}.ccac(comp),'og')
        
        title(["Eigenvalues/CanCorr Coeffs for both methods. Regularized: Green."])
    end
    % clear vars
    clear vars AUX d d0 d_long d0_long d_short d0_short t s REG_trn ADD_trn
    
end







