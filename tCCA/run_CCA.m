%% EXAMPLE SCRIPT
% BE AWARE:
% filtering and zscoring might be necessary/helpful. CCA components
% (regressors) have arbitrary scaling. Should not be a problem in the GLM
% regression, but keep it in mind!
clear all;
% ##### FOLLOWIG TWO LINES NEED CHANGE ACCRODING TO USER!
path.code = 'C:\Users\mayucel\Documents\PROJECTS\CODES\tCCA-GLM'; addpath(genpath(path.code)); % code directory
path.dir = 'C:\Users\mayucel\Google Drive\tCCA_GLM_PAPER\FB_RESTING_DATA'; % data directory
% #####
filename = 'resting_sim';
set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig))
set(groot,'defaultAxesCreateFcn',@(ax,~)set(ax.Toolbar,'Visible','off'))
sbjfolder = {'Subj33','Subj34','Subj36','Subj37','Subj38','Subj39', 'Subj40', 'Subj41', 'Subj43', 'Subj44','Subj46','Subj47','Subj49','Subj51'};


% options

flag_resting = 0; % will run on resting
flag_resting_half = 0; % use only the first half of resting to get mixing matrix.
rhoSD_ssThresh = 15;  % mm
stim_off = 1;  % for now no stim marks...
flag_plot = 1;
flag_save = 0;
HRFmin = -2;
HRFmax = 15; % used only for block design runs
Hb = 1; % 1 HbO / 2 HbR (for block only)
pre = 5;  % HRF range in sec to calculate ttest
post = 10;
segment = 20; % sec % segment for pwelch
% CCA parameters
param.tau = 2;% 2 (time points); %stepwidth for embedding in samples (tune to sample frequency!)
% param.NumOfEmb = 50 ; %75 200 12/ number of total embeddings (total time window will be tau * NumOfEmb)
param.ct = 0.5;   % correlation threshold
flags.pcaf =  [0 0]; % no pca of X or AUX
% timelag = 10%Y not full rank after ~ 30 (i.e. time shifted regressor number is higher than the number of time points) % in sec



tic;



for tl = 5% :1:10% time lag in sec
    
    timelag = tl
    
    for ss = 1:numel(sbjfolder) % loop across subjects
        ss
        cd([path.dir filesep sbjfolder{ss} filesep]);
        
        
        % load resting data for training
        [fq, t, AUX, d_long, d_short, d0_long, d0_short, d, d0, SD, s, lstLongAct,lstShortAct,lstHrfAdd] = load_nirs(filename);
        % AUX order : ACC1, ACC2, ACC3, PPG, BP, RESP
        param.NumOfEmb = timelag/param.tau * fq ;
        
        % merge data to right format
        AUX = [AUX, d0_short]; % full AUX = [acc1 acc2 acc3 PPG BP RESP, d_short];
        
        % check if the number of time points is odd/even, if odd make it even... (number of embedded should be the same)
        if mod(size(AUX,1),2) == 1
            AUX = AUX(1:end-1,:);
            d = d(1:end-1,:);
            d_long = d_long(1:end-1,:);
            d_short = d_short(1:end-1,:);
            d0 = d0(1:end-1,:);
            d0_long = d0_long(1:end-1,:);
            d0_short = d0_short(1:end-1,:);
            t = t(1:end-1,:);
            s = s(1:end-1,:);
        end
        
        
                % get first and second half
                c_half = round(size(AUX,1)/2);
        
            AUX1 = AUX(1:c_half,:);
            AUX2 = AUX(c_half+1:end,:);
            d_1 = d(1:c_half,:);
            d_2 = d(c_half+1:end,:);
            d_long_1 = d_long(1:c_half,:);
            d_long_2 = d_long(c_half+1:end,:);
            d_short_1 = d_short(1:c_half,:);
            d_short_2 = d_short(c_half+1:end,:);
            d0_1 = d0(1:c_half,:);
            d0_2 = d0(c_half+1:end,:);
            d0_long_1 = d0_long(1:c_half,:);
            d0_long_2 = d0_long(c_half+1:end,:);
            d0_short_1 = d0_short(1:c_half,:);
            d0_short_2 = d0_short(c_half+1:end,:);
            t_1 = t(1:c_half,:);
            t_2 = t(c_half+1:end,:);
            s_1 = s(1:c_half,:);
            s_2 = s(c_half+1:end,:);
        
        

    
        
        
        mean_AUX = mean(AUX1,1);
        AUX1 = AUX1 - repmat(mean_AUX, size(mean_AUX,1),1);
        mean_AUX = mean(AUX2,1);
        AUX2 = AUX2 - repmat(mean_AUX, size(mean_AUX,1),1);
        
        mean_d0_long_1 = mean(d0_long_1,1);
        X1 = d0_long_1 - repmat(mean_d0_long_1, size(d0_long_1,1),1);
        
        mean_d0_long_2 = mean(d0_long_2,1);
        X2 = d0_long_2 - repmat(mean_d0_long_2, size(d0_long_2,1),1);
        
        
        %% Perform the shiny script  % AUX = [acc1 acc2 acc3 PPG BP RESP, d_short];
        % Get ADD.Av from the training half
        [REG1,  ADD1] = perf_temp_emb_cca(X1,AUX1,param,flags);   % if first half training
        [REG2,  ADD2] = perf_temp_emb_cca(X2,AUX2,param,flags);  % if second half training
        
        % Get AUX (embedded) from the test half
        %% Temporally embed auxiliary data
        aux_sigs = AUX1;
        aux_emb1 = aux_sigs;
        for i=1:param.NumOfEmb
            aux=circshift( aux_sigs, i*param.tau, 1);
            aux(1:2*i,:)=repmat(aux(2*i+1,:),2*i,1);
            aux_emb1=[aux_emb1 aux];
        end
%         for i=1:param.NumOfEmb
%             aux=circshift( aux_sigs, -i*param.tau, 1);
%             aux(1:2*i,:)=repmat(aux(2*i+1,:),2*i,1);
%             aux_emb1=[aux_emb1 aux];
%         end
        
        
        
        aux_sigs = AUX2;
        aux_emb2 = aux_sigs;
        for i=1:param.NumOfEmb
            aux=circshift( aux_sigs, i*param.tau, 1);
            aux(1:2*i,:)=repmat(aux(2*i+1,:),2*i,1);
            aux_emb2=[aux_emb2 aux];
        end
%         for i=1:param.NumOfEmb
%             aux=circshift( aux_sigs, -i*param.tau, 1);
%             aux(1:2*i,:)=repmat(aux(2*i+1,:),2*i,1);
%             aux_emb2=[aux_emb2 aux];
%         end
        
        
        % Get Regressors
        V = aux_emb2;
        REG1 = (V - mean(V,1)) * ADD1.Av;   % first half training, second half test
        
        V = aux_emb1;
        REG2 = (V - mean(V,1)) * ADD2.Av;   % first half test, second half training
        
        
        % % take only the first ten regressors
        if size(REG1,2)>10
            Aaux1 = REG1(:,1:10);
        else
            Aaux1 = REG1;
        end
        
        if size(REG2,2)>10
            Aaux2 = REG2(:,1:10);
        else
            Aaux2 = REG2;
        end
        
        
        %% Add Regressors with high correlation to GLM
        window = ones(1,fq*segment); %rectangular window
        % window = hamming(segment*fs); % hamming
        xlim1 = 0;
        xlim2 = 25;
        
        dod = hmrIntensity2OD(d_2);
        dod = hmrBandpassFilt(dod, fq, 0, 0.5);
        dc_2 = hmrOD2Conc( dod, SD, [6 6]);
        
        
        dod = hmrIntensity2OD(d_1);
        dod = hmrBandpassFilt(dod, fq, 0, 0.5);
        dc_1 = hmrOD2Conc( dod, SD, [6 6]);
        
        %FIRST HALF TRAINING, SECOND HALF TESTING
        % GLM with SS
        [yavg_null, yavgstd, tHRF, nTrials, d_null, yresid, ysum2, beta, yR] = hmrDeconvHRF_DriftSS(dc_2, s_2, t_2, SD, [], [], [HRFmin HRFmax], 1, 1, [0.5 0.5], rhoSD_ssThresh, 1, 3, 0);

        % GLM with CCA outpout
        [yavg_new, yavgstd, tHRF, nTrials, d_new, yresid, ysum2, beta, yR] = hmrDeconvHRF_DriftSS(dc_2, s_2, t_2, SD, Aaux2, [], [HRFmin HRFmax], 1, 1, [0.5 0.5], 0, 0, 3, 0);
%         
        
%         %SECOND HALF TRAINING, FIRST HALF TESTING
%         % GLM with SS
%         [yavg_null, yavgstd, tHRF, nTrials, d_null, yresid, ysum2, beta, yR] = hmrDeconvHRF_DriftSS(dc_1, s_1, t_1, SD, [], [], [HRFmin HRFmax], 1, 1, [0.5 0.5], rhoSD_ssThresh, 1, 3, 0);
% 
%         % GLM with CCA outpout
%         [yavg_new, yavgstd, tHRF, nTrials, d_new, yresid, ysum2, beta, yR] = hmrDeconvHRF_DriftSS(dc_1, s_1, t_1, SD, Aaux1, [], [HRFmin HRFmax], 1, 1, [0.5 0.5], 0, 0, 3, 0);
% 
        
        
        
        
        lst_stim = find(s_1==1);
        lst_stim = lst_stim(1:nTrials);
        if lst_stim(1) < abs(HRFmin) * fq
            lst_stim = lst_stim(2:end);
        end
            
        for i = 1:size(lst_stim,1) % across trials
            HbO_null(:,i,:) = squeeze(d_null([lst_stim(i) - abs(HRFmin) * fq]:[lst_stim(i) + HRFmax * fq],Hb,:)); % get each trial (HbO)
            HbO_new(:,i,:) = squeeze(d_new([lst_stim(i) - abs(HRFmin) * fq]:[lst_stim(i) + HRFmax * fq],Hb,:)); % get each trial (HbO)
        end
        
        for i=1:size(HbO_null,3)
            MEAN_null(:,i)= nanmean(squeeze(HbO_null(:,:,i)),2);
            STD_null(:,i)=nanstd(squeeze(HbO_null(:,:,i)),0,2);
            MEAN_new(:,i)= nanmean(squeeze(HbO_new(:,:,i)),2);
            STD_new(:,i)=nanstd(squeeze(HbO_new(:,:,i)),0,2);
        end
        
        
        for i=1:size(HbO_null,3)
            MEAN_HRF_null(i,:)= nanmean(squeeze(HbO_null(pre*fq+abs(HRFmin*fq):post*fq+abs(HRFmin*fq),:,i)));  % channels by trials
            STD_HRF_null(i,:)=nanstd(squeeze(HbO_null(abs(HRFmin*fq):HRFmax*fq+abs(HRFmin*fq),:,i)));
            MEAN_baseline_null(i,:)= nanmean(squeeze(HbO_null(1:abs(HRFmin*fq),:,i)));
            STD_baseline_null(i,:)=nanstd(squeeze(HbO_null(1:abs(HRFmin*fq),:,i)));
            
            MEAN_HRF_new(i,:)= nanmean(squeeze(HbO_new(pre*fq+abs(HRFmin*fq):post*fq+abs(HRFmin*fq),:,i)));  % channels by trials
            STD_HRF_new(i,:)=nanstd(squeeze(HbO_new(pre*fq+abs(HRFmin*fq):post*fq+abs(HRFmin*fq),:,i)));
            MEAN_baseline_new(i,:)= nanmean(squeeze(HbO_new(1:abs(HRFmin*fq),:,i)));
            STD_baseline_new(i,:)=nanstd(squeeze(HbO_new(1:abs(HRFmin*fq),:,i)));
        end
        
        for i=1:size(HbO_null,3)
            [h,p,c,stats]=ttest(MEAN_HRF_null(i,:),(MEAN_baseline_null(i,:)));
            pOxy_null(ss,i)=p;
            [h,p,c,stats]=ttest(MEAN_HRF_new(i,:),(MEAN_baseline_new(i,:)));
            pOxy_new(ss,i)=p;
        end
        
        p_null = zeros(size(pOxy_null));
        p_new = zeros(size(pOxy_new));
        
        % get number of active channels
        for i = 1:size(pOxy_new,1)  % indice includes both long and short
            p_null(i,find(pOxy_null(i,:)<=0.05)) = 1;
            p_new(i,find(pOxy_new(i,:)<=0.05)) = 1;
        end
        
        ml = SD.MeasList;
        mlAct = SD.MeasListAct;
        lst = find(ml(:,4)==1);
        rhoSD = zeros(length(lst),1);
        posM = zeros(length(lst),3);
        for iML = 1:length(lst)
            rhoSD(iML) = sum((SD.SrcPos(ml(lst(iML),1),:) - SD.DetPos(ml(lst(iML),2),:)).^2).^0.5;
            posM(iML,:) = (SD.SrcPos(ml(lst(iML),1),:) + SD.DetPos(ml(lst(iML),2),:)) / 2;
        end
        lstSS = lst(find(rhoSD<15));
        p_null(:,lstSS)= 0;  % remove ss channels from list
        p_new(:,lstSS)= 0;
        
        
        if flag_plot
            plot_block(MEAN_null, MEAN_new, HRFmin, HRFmax, fq, pOxy_null, pOxy_new, ss, STD_null,STD_new, tHRF,timelag,path.dir,lstHrfAdd);
        end
        
        clear MEAN_HRF_null MEAN_HRF_new   MEAN_baseline_null  MEAN_baseline_new STD_HRF_null STD_HRF_new STD_baseline_null STD_baseline_new...
            HbO_null HbO_new  MEAN_null  MEAN_new  STD_null  STD_new REG ADD X Aaux d dod dc s t AUX lst_stim h p c stats
        
        
        
    end
end


toc;

