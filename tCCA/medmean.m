function [C,M,P,F] = medmean(CORR, MSE, PVAL, FSCORE, mflag)
%MEDMEAN Summary of this function goes here

% Validation parameters
tlags = 0:1:10;
stpsize = 2:2:24;
cthresh = 0:0.1:0.9;


switch mflag
%% MEAN
    case 1
        %% Average across channels
        C = squeeze(nanmean(CORR,2));
        M = squeeze(median(MSE,2));
        P = squeeze(median(PVAL,2));
        %% now average across splits
        C = squeeze(median(C,3));
        M = squeeze(median(M,3));
        P = squeeze(median(P,3));
        F = squeeze(median(FSCORE,3));
        %% now average across subjects
        C = median(C,1);
        M = median(M,1);
        P = median(P,1);
        F = median(F,1);
        %% MEDIAN
    case 2
        %% Append all subjects, channels and folds
        for tt = 1:numel(tlags)
            for ss = 1:numel(stpsize)
                for cc = 1:numel(cthresh)
                    for hh = 1:2
                        buf = CORR(:,:,hh,:,tt,ss,cc);
                        C(:,hh,tt,ss,cc) = buf(:);
                        buf = MSE(:,:,hh,:,tt,ss,cc);
                        M(:,hh,tt,ss,cc) = buf(:);
                        buf = FSCORE(:,hh,:,tt,ss,cc);
                        F(:,hh,tt,ss,cc) = buf(:);
                        buf = PVAL(:,:,hh,:,tt,ss,cc);
                        P(:,hh,tt,ss,cc) = buf(:);
                    end
                end
            end
        end
        %% median across across channels, subjects and folds
        C = median(C,1,'omitnan');
        M = median(M,1,'omitnan');
        P = median(P,1,'omitnan');
        F = median(F,1,'omitnan');
end
end

