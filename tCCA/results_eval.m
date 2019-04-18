function [p_SS, p_CCA, pOxy_SS, pOxy_CCA] = results_eval(sbj, d_ss, d_cca, tHRF, timelag, lst_stim, SD, fq, lstHrfAdd, eval_param, flag_plot, path)
%PLOT_EVAL Summary of this function goes here


for i = 1:size(lst_stim,1) % across trials
    HbO_SS(:,i,:) = squeeze(d_ss([lst_stim(i) - abs(eval_param.HRFmin) * fq]:[lst_stim(i) + eval_param.HRFmax * fq],eval_param.Hb,:)); % get each trial (HbO)
    HbO_CCA(:,i,:) = squeeze(d_cca([lst_stim(i) - abs(eval_param.HRFmin) * fq]:[lst_stim(i) + eval_param.HRFmax * fq],eval_param.Hb,:)); % get each trial (HbO)
end

for i=1:size(HbO_SS,3)
    MEAN_SS(:,i)= nanmean(squeeze(HbO_SS(:,:,i)),2);
    STD_SS(:,i)=nanstd(squeeze(HbO_SS(:,:,i)),0,2);
    MEAN_CCA(:,i)= nanmean(squeeze(HbO_CCA(:,:,i)),2);
    STD_SS(:,i)=nanstd(squeeze(HbO_CCA(:,:,i)),0,2);
end


for i=1:size(HbO_SS,3)
    MEAN_HRF_SS(i,:)= nanmean(squeeze(HbO_SS(eval_param.pre*fq+abs(eval_param.HRFmin*fq):eval_param.post*fq+abs(eval_param.HRFmin*fq),:,i)));  % channels by trials
    STD_HRF_SS(i,:)=nanstd(squeeze(HbO_SS(abs(eval_param.HRFmin*fq):eval_param.HRFmax*fq+abs(eval_param.HRFmin*fq),:,i)));
    MEAN_baseline_SS(i,:)= nanmean(squeeze(HbO_SS(1:abs(eval_param.HRFmin*fq),:,i)));
    STD_baseline_SS(i,:)=nanstd(squeeze(HbO_SS(1:abs(eval_param.HRFmin*fq),:,i)));
    
    MEAN_HRF_CCA(i,:)= nanmean(squeeze(HbO_CCA(eval_param.pre*fq+abs(eval_param.HRFmin*fq):eval_param.post*fq+abs(eval_param.HRFmin*fq),:,i)));  % channels by trials
    STD_HRF_CCA(i,:)=nanstd(squeeze(HbO_CCA(eval_param.pre*fq+abs(eval_param.HRFmin*fq):eval_param.post*fq+abs(eval_param.HRFmin*fq),:,i)));
    MEAN_baseline_CCA(i,:)= nanmean(squeeze(HbO_CCA(1:abs(eval_param.HRFmin*fq),:,i)));
    STD_baseline_CCA(i,:)=nanstd(squeeze(HbO_CCA(1:abs(eval_param.HRFmin*fq),:,i)));
end

for i=1:size(HbO_SS,3)
    format long
    [h,p,c,stats]=ttest(MEAN_HRF_SS(i,:),(MEAN_baseline_SS(i,:)));
    pOxy_SS(sbj,i)=p;
    [h,p,c,stats]=ttest(MEAN_HRF_CCA(i,:),(MEAN_baseline_CCA(i,:)));
    pOxy_CCA(sbj,i)=p;
end

p_SS = zeros(size(pOxy_SS,2),1);
p_CCA = zeros(size(pOxy_CCA,2),1);

% get number of active channels that we added HRF (True Positives only!)
lst = find(pOxy_SS(sbj,:)<=0.05);
foo = ismember(lst,lstHrfAdd(:,1));
lst = lst(find(foo==1));
p_SS(lst) = 1;

lst = find(pOxy_CCA(sbj,:)<=0.05);
foo = ismember(lst,lstHrfAdd(:,1));
lst = lst(find(foo==1));
p_CCA(lst) = 1;

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
p_SS(lstSS)= 0;  % remove ss channels from list
p_CCA(lstSS)= 0;


if flag_plot
    plot_block(MEAN_SS, MEAN_CCA, eval_param.HRFmin, eval_param.HRFmax, fq, pOxy_SS, pOxy_CCA, sbj, STD_SS,STD_SS, tHRF, timelag, path.dir,lstHrfAdd);
end

% average pvals of activated channels
buf = pOxy_SS(sbj,:);
pvals_SS = buf(find(p_SS));
disp(['avg p SS ' num2str(nanmean(pvals_SS))])
buf = pOxy_CCA(sbj,:);
pvals_CCA = buf(find(p_CCA));
disp(['avg p CCA ' num2str(nanmean(pvals_CCA))])

%% lets do the calculation of the other performance metrics here later on
%
%
%

end

