function [pval_SS, pval_CCA] = results_eval_block(sbj, d_ss, d_cca, yavg_ss, yavg_cca, tHRF, timelag, sts, ctr, lst_stim, SD, fq, eval_param, flag_plot, path)
%PLOT_EVAL Summary of this function goes here

% get hrf at each trial
for j = 1:2 % across HbO/HbR
    for i = 1:size(lst_stim,1) % across trials
        Hb_SS(:,i,:,j) = squeeze(d_ss([lst_stim(i) - abs(eval_param.HRFmin) * fq]:[lst_stim(i) + eval_param.HRFmax * fq],j,:)); % get each trial (HbO)
        Hb_CCA(:,i,:,j) = squeeze(d_cca([lst_stim(i) - abs(eval_param.HRFmin) * fq]:[lst_stim(i) + eval_param.HRFmax * fq],j,:)); % get each trial (HbO)
    end
end
% Hb_SS: # of time points X # of trials X # of channels X HbO/HbR


% take the mean across between eval_param.pre and post generally around peak HRF
for j = 1:2 % HbO/HbR
    for i=1:size(Hb_SS,3) % across channels
        % HBO & HBR
        MEAN_HRF_SS(i,:,j)= nanmean(squeeze(Hb_SS(eval_param.pre*fq+abs(eval_param.HRFmin*fq):eval_param.post*fq+abs(eval_param.HRFmin*fq),:,i,j)));  % channels by trials
        MEAN_baseline_SS(i,:,j)= nanmean(squeeze(Hb_SS(1:abs(eval_param.HRFmin*fq),:,i,j)));
        
        MEAN_HRF_CCA(i,:,j)= nanmean(squeeze(Hb_CCA(eval_param.pre*fq+abs(eval_param.HRFmin*fq):eval_param.post*fq+abs(eval_param.HRFmin*fq),:,i,j)));  % channels by trials
        MEAN_baseline_CCA(i,:,j)= nanmean(squeeze(Hb_CCA(1:abs(eval_param.HRFmin*fq),:,i,j)));
    end
end

% get stats
for i=1:size(Hb_SS,3)
    format long
    %HbO & HbR
    for ii=1:2
        [h,p,c,stats]=ttest(MEAN_HRF_SS(i,:,ii),(MEAN_baseline_SS(i,:,ii)));
        pval_SS(i,ii)=p;
        [h,p,c,stats]=ttest(MEAN_HRF_CCA(i,:,ii),(MEAN_baseline_CCA(i,:,ii)));
        pval_CCA(i,ii)=p;
    end
end


% get mean across trials for the plot_block function
for i=1:size(Hb_SS,3)
    for j = 1:2
        % HbO and HbR
        MEAN_SS(:,i,j)= nanmean(squeeze(Hb_SS(:,:,i,j)),2);
        STD_SS(:,i,j)=nanstd(squeeze(Hb_SS(:,:,i,j)),0,2);
        MEAN_CCA(:,i,j)= nanmean(squeeze(Hb_CCA(:,:,i,j)),2);
        STD_CCA(:,i,j)=nanstd(squeeze(Hb_CCA(:,:,i,j)),0,2);
        
    end
end



if flag_plot
    plot_visual_block(MEAN_SS, MEAN_CCA, eval_param.HRFmin, eval_param.HRFmax, fq, pval_SS, pval_CCA, sbj, STD_SS,STD_SS, tHRF, timelag, sts, ctr, path.dir);
end



end

