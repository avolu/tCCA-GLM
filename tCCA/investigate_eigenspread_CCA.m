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
    path.code = 'E:\Office\Research\Software - Scripts\Matlab\Regression tCCA GLM\tCCA-GLM'; addpath(genpath(path.code)); % code directory
    path.dir = 'C:\Users\mladm\Google Drive\tCCA_GLM_PAPER\FB_RESTING_DATA'; % data directory
    path.save = 'C:\Users\mladm\Google Drive\tCCA_GLM_PAPER'; % save directory
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
        X = d0_long(trnIDX,:);
        % use the test data of aux
        Y = aux_emb;
        
        %% This is where the following function would run        
        [REG_trn{tt},  ADD_trn{tt}] = perf_temp_emb_cca(X,AUX(trnIDX,:),param,flags);
        %% But we formulate the generalised eigenvalue equation to investigate the eigenvalue spectrum
        
        % Center the variables
        X = X - repmat(mean(X,1), size(X,1), 1);
        Y = Y - repmat(mean(Y,1), size(X,1), 1);
        % time points
        T = size(X,1);
        % Calculate empirical Auto-Covariance matrices
        Cxx = X'*X/T;
        Cyy = Y'*Y/T;
        % Calculate empirical Cross-Covariance matrices
        Cxy = X'*Y/T;
        Cyx = Y'*X/T;
        
        % put auto-covariance matrices into block matrix form (generalized eigenvalue problem) 
        [dx, dy] = size(Cxy);
        % A = [0 Cxy; Cyx 0]
        A = [zeros(dx,dx) Cxy; Cyx zeros(dy, dy)];
        % B = [Cxx 0; 0 Cyy]
        B = [Cxx zeros(dx,dy); zeros(dy, dx) Cyy ];

        % add regularization 
        gamma = 0%.01;
        B = B + gamma*eye(dx + dy);
        
        % calculate eigenvalues of B
        lambda_b = eig(B);
        figure
        plot(lambda_b)
        
        % Calculate Generalized Eigenvalues
        %   [V,D] = EIG(A,B) produces a diagonal matrix D of generalized
        %   eigenvalues and a full matrix V whose columns are the corresponding
        %   eigenvectors so that A*V = B*V*D.
        [V, D] = eig(A,B);
        lambda=diag(D);
        figure
        plot(lambda)
        
        Wx = V(1:dx, 1:dx);
        Wy = V(dx+1:end, dx+1:end);
       
        Sx = X*Wx;
        Sy = Y*Wy;
        
    end
% clear vars
clear vars AUX d d0 d_long d0_long d_short d0_short t s REG_trn ADD_trn

end







