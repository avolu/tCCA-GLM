function [DET_SS, DET_CCA, pval_SS, pval_CCA, ROCLAB, MSE_SS, MSE_CCA, CORR_SS, CORR_CCA] = results_eval(sbj, d_ss, d_cca, tHRF, timelag, lst_stim, SD, fq, lstHrfAdd, eval_param, flag_plot, path, hrf)
%PLOT_EVAL Summary of this function goes here

% find short separation channels
ml = SD.MeasList;
mlAct = SD.MeasListAct;
lst_sig = find(ml(:,4)==1);
rhoSD = zeros(length(lst_sig),1);
posM = zeros(length(lst_sig),3);
for iML = 1:length(lst_sig)
    rhoSD(iML) = sum((SD.SrcPos(ml(lst_sig(iML),1),:) - SD.DetPos(ml(lst_sig(iML),2),:)).^2).^0.5;
    posM(iML,:) = (SD.SrcPos(ml(lst_sig(iML),1),:) + SD.DetPos(ml(lst_sig(iML),2),:)) / 2;
end
lstSS = lst_sig(find(rhoSD<15));

for i = 1:size(lst_stim,1) % across trials
    HbO_SS(:,i,:) = squeeze(d_ss([lst_stim(i) - abs(eval_param.HRFmin) * fq]:[lst_stim(i) + eval_param.HRFmax * fq],1,:)); % get each trial (HbO)
    HbO_CCA(:,i,:) = squeeze(d_cca([lst_stim(i) - abs(eval_param.HRFmin) * fq]:[lst_stim(i) + eval_param.HRFmax * fq],1,:)); % get each trial (HbO)
    HbR_SS(:,i,:) = squeeze(d_ss([lst_stim(i) - abs(eval_param.HRFmin) * fq]:[lst_stim(i) + eval_param.HRFmax * fq],2,:)); % get each trial (HbR)
    HbR_CCA(:,i,:) = squeeze(d_cca([lst_stim(i) - abs(eval_param.HRFmin) * fq]:[lst_stim(i) + eval_param.HRFmax * fq],2,:)); % get each trial (HbR)
end

for i=1:size(HbO_SS,3)
    % HBO
    MEAN_SS(:,i,1)= nanmean(squeeze(HbO_SS(:,:,i)),2);
    STD_SS(:,i,1)=nanstd(squeeze(HbO_SS(:,:,i)),0,2);
    MEAN_CCA(:,i,1)= nanmean(squeeze(HbO_CCA(:,:,i)),2);
    STD_SS(:,i,1)=nanstd(squeeze(HbO_CCA(:,:,i)),0,2);
    % HBR
    MEAN_SS(:,i,2)= nanmean(squeeze(HbR_SS(:,:,i)),2);
    STD_SS(:,i,2)=nanstd(squeeze(HbR_SS(:,:,i)),0,2);
    MEAN_CCA(:,i,2)= nanmean(squeeze(HbR_CCA(:,:,i)),2);
    STD_SS(:,i,2)=nanstd(squeeze(HbR_CCA(:,:,i)),0,2);
end


for i=1:size(HbO_SS,3) % across channels
    % HBO & HBR
    MEAN_HRF_SS(i,:,1)= nanmean(squeeze(HbO_SS(eval_param.pre*fq+abs(eval_param.HRFmin*fq):eval_param.post*fq+abs(eval_param.HRFmin*fq),:,i)));  % channels by trials
    MEAN_baseline_SS(i,:,1)= nanmean(squeeze(HbO_SS(1:abs(eval_param.HRFmin*fq),:,i)));
    
    MEAN_HRF_CCA(i,:,1)= nanmean(squeeze(HbO_CCA(eval_param.pre*fq+abs(eval_param.HRFmin*fq):eval_param.post*fq+abs(eval_param.HRFmin*fq),:,i)));  % channels by trials
    MEAN_baseline_CCA(i,:,1)= nanmean(squeeze(HbO_CCA(1:abs(eval_param.HRFmin*fq),:,i)));
    % HBR
    MEAN_HRF_SS(i,:,2)= nanmean(squeeze(HbR_SS(eval_param.pre*fq+abs(eval_param.HRFmin*fq):eval_param.post*fq+abs(eval_param.HRFmin*fq),:,i)));  % channels by trials
    MEAN_baseline_SS(i,:,2)= nanmean(squeeze(HbR_SS(1:abs(eval_param.HRFmin*fq),:,i)));
    
    MEAN_HRF_CCA(i,:,2)= nanmean(squeeze(HbR_CCA(eval_param.pre*fq+abs(eval_param.HRFmin*fq):eval_param.post*fq+abs(eval_param.HRFmin*fq),:,i)));  % channels by trials
    MEAN_baseline_CCA(i,:,2)= nanmean(squeeze(HbR_CCA(1:abs(eval_param.HRFmin*fq),:,i)));
end

for i=1:size(HbO_SS,3)
    format long
    %HbO & HbR
    for ii=1:2
        [h,p,c,stats]=ttest(MEAN_HRF_SS(i,:,ii),(MEAN_baseline_SS(i,:,ii)));
        pval_SS(i,ii)=p;
        [h,p,c,stats]=ttest(MEAN_HRF_CCA(i,:,ii),(MEAN_baseline_CCA(i,:,ii)));
        pval_CCA(i,ii)=p;
    end
end

DET_SS = zeros(size(MEAN_SS,2),2);
DET_CCA = zeros(size(MEAN_CCA,2),2);

% get number of active channels that we added HRF (True Positives only!)
% list of true negatives
nohrflist = 1:1:size(HbO_SS,3);
rmvidx = [lstHrfAdd(:,1); lstSS];
nohrflist(rmvidx)=[];
% init variables
TP_SS=NaN(16,2);
FP_SS=NaN(16,2);
TN_SS=NaN(18,2);
TN_CCA=NaN(18,2);
TP_CCA=NaN(16,2);
FP_CCA=NaN(16,2);

%HbO & HbR
for ii=1:2
    % SS
    lst_sig = find(pval_SS(:,ii)<=0.05);
    TP_SS(1:numel(lst_sig),ii) = ismember(lst_sig,lstHrfAdd(:,1)); % # of true positive channels
    FP_SS(1:numel(lst_sig),ii) = ~ismember(lst_sig,lstHrfAdd(:,1)); % # of false positive channels
    DET_SS(lst_sig(find(TP_SS(:,ii)==1)),ii) = 1;  % TP
    DET_SS(lst_sig(find(FP_SS(:,ii)==1)),ii) = -1; % FP
    lst_notsig = find(pval_SS(:,ii)>0.05);
    TN_SS(1:numel(lst_notsig),ii) = ismember(lst_notsig,nohrflist);
    FN_SS(1:numel(lst_notsig),ii) = ismember(lst_notsig,lstHrfAdd(:,1));
    DET_SS(lst_notsig(find(FN_SS(:,ii)==1)),ii) = 2; % FN
    DET_SS(lst_notsig(find(TN_SS(:,ii)==1)),ii) = -2; % TN
    
    % CCA
    lst_sig = find(pval_CCA(:,ii)<=0.05);
    TP_CCA(1:numel(lst_sig),ii) = ismember(lst_sig,lstHrfAdd(:,1)); % # of true positive channels
    FP_CCA(1:numel(lst_sig),ii) = ~ismember(lst_sig,lstHrfAdd(:,1)); % # of false positive channels
    DET_CCA(lst_sig(find(TP_CCA(:,ii)==1)),ii) = 1; % TP
    DET_CCA(lst_sig(find(FP_CCA(:,ii)==1)),ii) = -1; % FP
    lst_notsig = find(pval_CCA(:,ii)>0.05);
    TN_CCA(1:numel(lst_notsig),ii) = ismember(lst_notsig,nohrflist);
    FN_CCA(1:numel(lst_notsig),ii) = ismember(lst_notsig,lstHrfAdd(:,1));
    DET_CCA(lst_notsig(find(FN_CCA(:,ii)==1)),ii) = 2; % FN
    DET_CCA(lst_notsig(find(TN_CCA(:,ii)==1)),ii) = -2; % TN
    
    ROCLAB.val = [1,-1,2,-2,0];
    ROCLAB.name = {'TP','FP','FN','TN', 'PRND'};
end

% remove ss channels from list
DET_SS(lstSS,:)= 0;
DET_CCA(lstSS,:)= 0;


%% lets do the calculation of the other performance metrics here later on
% correlation
% cut to the same timebase and baseline correct
MEAN_SS_ev = MEAN_SS(abs(eval_param.HRFmin*fq)-1:end,:,:);
MEAN_SS_ev = MEAN_SS_ev-nanmean(MEAN_SS(1:abs(fq*eval_param.HRFmin),:,:),1);
MEAN_CCA_ev = MEAN_CCA(abs(eval_param.HRFmin*fq)-1:end,:,:);
MEAN_CCA_ev = MEAN_CCA_ev-nanmean(MEAN_CCA(1:abs(fq*eval_param.HRFmin),:,:),1);
hrfeval = hrf.hrf_conc(1:size(MEAN_SS_ev,1),:);

%init variables
CORR_SS=NaN(16,2);
CORR_CCA=NaN(16,2);
MSE_SS=NaN(16,2);
MSE_CCA=NaN(16,2);

for ii=1:2
    idx = lstHrfAdd(:,1);
    
    % calculate correlation and MSE
    buf = corr(squeeze(MEAN_SS_ev(:,idx,ii)),squeeze(hrfeval(:,ii)));
    CORR_SS(1:numel(buf),ii) = buf;
    buf = corr(squeeze(MEAN_CCA_ev(:,idx,ii)),squeeze(hrfeval(:,ii)));
    CORR_CCA(1:numel(buf),ii) = buf;
    buf = sqrt(nanmean((squeeze(MEAN_SS_ev(:,idx,ii))-squeeze(hrfeval(:,ii))).^2));
    MSE_SS(1:numel(buf),ii) = buf;
    buf = sqrt(nanmean((squeeze(MEAN_CCA_ev(:,idx,ii))-squeeze(hrfeval(:,ii))).^2));
    MSE_CCA(1:numel(buf),ii) = buf;
end


if flag_plot
    plot_block(MEAN_SS, MEAN_CCA, CORR_SS, CORR_CCA, MSE_SS, MSE_CCA, eval_param.HRFmin, eval_param.HRFmax, fq, pval_SS, pval_CCA, sbj, STD_SS,STD_SS, tHRF, timelag, path.dir,lstHrfAdd,hrf);
end



end

