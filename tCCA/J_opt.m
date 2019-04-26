function [fval] = J_opt(CORR, MSE, PVAL, FSCORE, Jparam, reg)
% J_OPT Cost function to optimize
% uses results from run_CCA
% CORR, MSE, PVAL, FSCORE, CHNO are 4D matrices with results from run_CCA
% averaged across subjects, folds and channels:
%   (HbO/R (2) x tlags (11) x stpsize (12) x cthresh (10))
%
%% INPUTS

% Jparam: struct with parameters/options for optimization
%   .fact:      struct with factors (weights) for adapting J
%       .corr
%       .mse
%       .pval
%       .fscore
%       .HbO
%       .HbR
%   .mtype:     optimization using (subjects x channels x folds)
%        = 1 mean
%        = 2 all
%        = 3 median
%   .nflag:     normalization approach
%       1: X/max
%       2:(X-min)/(max-min)
%   .thresh:    segmentation approach: threshold for segmentation
% reg: take neighboring points in all directions into account
%       = 0 (off) or any other number (on)

% Validation parameters
tlags = 0:1:10;
stpsize = 2:2:24;
cthresh = 0:0.1:0.9;


%% MEDIAN or MEAN or NOTHING across channels, subjects and folds
% #SBJ x #CH x 2(Hbo+HbR) x 2 (cv split) x tlag x stepsize x corrthres
switch Jparam.mtype
    %% MEAN
    case 1
        [C,M,P,F] = medmean(CORR, MSE, PVAL, FSCORE, 1);
        %% MEDIAN
    case 2
        [C,M,P,F] = medmean(CORR, MSE, PVAL, FSCORE, 2);
        %% ALL
    case 3
        [C,M,P,F] = medmean(CORR, MSE, PVAL, FSCORE, 3);
end


%% NORMALIZE METRICS FOR OBJECTIVE FUNCTION
switch Jparam.nflag
    %% X/(max(X)
    case 1
        for hh=1:2
            C(:,hh,:,:,:) = C(:,hh,:,:,:)./max(C(:,hh,:));
            M(:,hh,:,:,:) = M(:,hh,:,:,:)./max(M(:,hh,:));
            P(:,hh,:,:,:) = P(:,hh,:,:,:)./max(P(:,hh,:));
            F(:,hh,:,:,:) = F(:,hh,:,:,:)./max(F(:,hh,:));
        end
        %% (X-min(X))/(max(X)-min(X))
    case 2
        for hh=1:2
            for pp = 1:size(C,1)
                C(pp,hh,:,:,:) = (C(pp,hh,:,:,:)-min(squeeze(C(pp,hh,:))))/(max(squeeze(C(pp,hh,:)))-min(squeeze(C(pp,hh,:))));
                M(pp,hh,:,:,:) = (M(pp,hh,:,:,:)-min(squeeze(M(pp,hh,:))))/(max(squeeze(M(pp,hh,:)))-min(squeeze(M(pp,hh,:))));
                P(pp,hh,:,:,:) = (P(pp,hh,:,:,:)-min(squeeze(P(pp,hh,:))))/(max(squeeze(P(pp,hh,:)))-min(squeeze(P(pp,hh,:))));
            end
            for pp = 1:size(F,1)
                F(pp,hh,:,:,:) = (F(pp,hh,:,:,:)-min(squeeze(F(pp,hh,:))))/(max(squeeze(F(pp,hh,:)))-min(squeeze(F(pp,hh,:))));
            end
        end
end


%% Calculate Objective Function Output
for tt = 1:numel(tlags)
    for ss = 1:numel(stpsize)
        for cc = 1:numel(cthresh)
            fval(tt,ss,cc) = 0;
            % get region indices, throw points outside of parameter set
            tr = [tt-reg.step:tt+reg.step];
            tr(find(tr<1))=[];
            tr(find(tr>numel(tlags)))=[];
            sr = [ss-reg.step:ss+reg.step];
            sr(find(sr<1))=[];
            sr(find(sr>numel(stpsize)))=[];
            cr = [cc-reg.step:cc+reg.step];
            cr(find(cr<1))=[];
            cr(find(cr>numel(cthresh)))=[];
            % determine divider (number of points in region without origin)
            div = numel(tr)*numel(sr)*numel(cr)-1;
            % calc at points and weighted region around that point
            for ti = tr
                for si = sr
                    for ci = cr
                        if ti == tt && si == ss && ci == cc
                            divider = 1;
                        else
                            divider = div/reg.weight;
                        end
                        fval(tt,ss,cc) = fval(tt,ss,cc) + ...
                            (- sum(Jparam.fact.HbO*Jparam.fact.corr*C(:,1,ti,si,ci), 'omitnan') ...
                            - sum(Jparam.fact.HbR*Jparam.fact.corr*C(:,2,ti,si,ci), 'omitnan') ...
                            + sum(Jparam.fact.HbO*Jparam.fact.mse*M(:,1,ti,si,ci), 'omitnan') ...
                            + sum(Jparam.fact.HbR*Jparam.fact.mse*M(:,2,ti,si,ci), 'omitnan') ...
                            + sum(Jparam.fact.HbO*Jparam.fact.pval*P(:,1,ti,si,ci), 'omitnan') ...
                            + sum(Jparam.fact.HbR*Jparam.fact.pval*P(:,2,ti,si,ci), 'omitnan') ...
                            - sum(Jparam.fact.HbO*Jparam.fact.fscore*F(:,1,ti,si,ci), 'omitnan') ...
                            - sum(Jparam.fact.HbR*Jparam.fact.fscore*F(:,2,ti,si,ci), 'omitnan')) ...
                            /divider;
                    end
                end
            end
        end
    end
end



%% only for median/mean approach
if Jparam.mtype == 1 || Jparam.mtype == 2
    %% Global optimum using segmentation: find parameter set(s) that are in each optimal segment
    % and find overlaps
    % HbO & HbR
    for hh=1:2
        buf = squeeze(C(:,hh,:,:,:));
        [t,s,c] = ind2sub(size(buf),find(buf>Jparam.thresh));
        Cs = [t s c];
        buf = squeeze(M(:,hh,:,:,:));
        [t,s,c] = ind2sub(size(buf),find(buf<(1-Jparam.thresh)));
        Ms = [t s c];
        buf = squeeze(F(:,hh,:,:,:));
        [t,s,c] = ind2sub(size(buf),find(buf>Jparam.thresh));
        Fs = [t s c];
        
        R1 =intersect(Cs,Ms,'rows');
        R{hh} = intersect(R1,Fs,'rows');
    end
    figure
    for ii=1:10
        subplot(3,4,ii)
        idx = find(R{hh}(:,3) == ii);
        plot(stpsize(R{1}(idx,2)),tlags(R{1}(idx,1)), 'or')
        hold on
        plot(stpsize(R{2}(idx,2)),tlags(R{2}(idx,1)), '+b')
        xlabel('stepsize / smpl')
        ylabel('time lags / s')
        xlim([2 24])
        ylim([0 10])
        title(['optimal intersecting parameters HbO/HbR, ctrsh: ' num2str(cthresh(ii))])
    end
end

end

