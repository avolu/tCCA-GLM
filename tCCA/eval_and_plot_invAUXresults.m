clear all

%% +++++++++++++++++++++++
%% SCRIPT CONFIGURATION
% +++++++++++++++++++++++
% user: 1 Meryem | 0 Alex
melexflag = 1;
% select which hrf amplitude data: 1 (20%), 2 (50%) or 3 (100%)
hhh = [2];
% select which metric type: 1 (average of single trial HRF MSEs), 2: MSE of average HRF
mmm = [1];
% Use only true positives for evaluation of metrics
TP_flag = true;
% indices of optimal parameterset
pOptfix = [4 7 6];



%% Data
% ##### FOLLOWING TWO LINES NEED CHANGE ACCORDING TO USER!
if melexflag
    %Meryem
    path.code = 'C:\Users\mayucel\Documents\PROJECTS\CODES\tCCA-GLM'; addpath(genpath(path.code)); % code directory
    path.dir = 'C:\Users\mayucel\Google Drive\tCCA_GLM_PAPER\FB_RESTING_DATA'; % data directory
    path.save = 'C:\Users\mayucel\Google Drive\tCCA_GLM_PAPER'; % save directory
    %path.auxres20 = 'C:\Users\mayucel\Google Drive\tCCA_GLM_PAPER\CV_AUX_contributions_20'; % save directory
    %path.auxres50 = 'C:\Users\mayucel\Google Drive\tCCA_GLM_PAPER\CV_AUX_contributions_50'; % save directory
    %path.auxres100 = 'C:\Users\mayucel\Google Drive\tCCA_GLM_PAPER\CV_AUX_contributions_100'; % save directory
    %path.auxres20stmse = 'C:\Users\mayucel\Google Drive\tCCA_GLM_PAPER\CV_AUX_contributions_20_stMSE'; % save directory
    path.auxres50stmse = 'C:\Users\mayucel\Google Drive\tCCA_GLM_PAPER\CV_AUX_contributions_50_stMSE'; % save directory
    %path.auxres100stmse = 'C:\Users\mayucel\Google Drive\tCCA_GLM_PAPER\CV_AUX_contributions_100_stMSE'; % save directory
else
    %Alex
    path.code = 'D:\Office\Research\Software - Scripts\Matlab\Regression tCCA GLM\tCCA-GLM'; addpath(genpath(path.code)); % code directory
    path.dir = 'C:\Users\avolu\Google Drive\tCCA_GLM_PAPER\FB_RESTING_DATA'; % data directory
    path.save = 'C:\Users\avolu\Google Drive\tCCA_GLM_PAPER'; % save directory
    %path.auxres20 = 'C:\Users\avolu\Google Drive\tCCA_GLM_PAPER\CV_AUX_contributions_20'; % save directory
    %path.auxres50 = 'C:\Users\avolu\Google Drive\tCCA_GLM_PAPER\CV_AUX_contributions_50'; % save directory
    %path.auxres100 = 'C:\Users\avolu\Google Drive\tCCA_GLM_PAPER\CV_AUX_contributions_100'; % save directory
    %path.auxres20stmse = 'C:\Users\avolu\Google Drive\tCCA_GLM_PAPER\CV_AUX_contributions_20_stMSE'; % save directory
    path.auxres50stmse = 'C:\Users\avolu\Google Drive\tCCA_GLM_PAPER\CV_AUX_contributions_50_stMSE'; % save directory
    %path.auxres100stmse = 'C:\Users\avolu\Google Drive\tCCA_GLM_PAPER\CV_AUX_contributions_100_stMSE'; % save directory
end

% Validation parameters
tlags = 0:1:10;
stpsize = 2:2:24;
cthresh = 0:0.1:0.9;
evparams.tlags = tlags;
evparams.stpsize = stpsize;
evparams.cthresh = cthresh;

hblab = {'HbO', 'HbR'};
metrttl = {'single trial', 'block avg'};

for metr=mmm
    for hrff=hhh
        %% load results data from all subjects
        % Dimensions of output metrics
        % # of sbjs x #CH x 2(Hbo+HbR) x 2 (cv split) x tlag x stepsize x corrthres
            switch hrff
                case 1
                    hrfamp = 20;
                    switch metr
                        case 1
                            load([path.auxres20stmse  '\results_all_sbj.mat']);
                        case 2
                            load([path.auxres20 '\results_all_sbj.mat']);         
                    end
                case 2
                    hrfamp = 50;
                    switch metr
                        case 1
                            load([path.auxres50stmse  '\results_all_sbj.mat']);
                        case 2
                            load([path.auxres50 '\results_all_sbj.mat']);
                    end
                case 3
                    hrfamp = 100;
                    switch metr
                        case 1
                            load([path.auxres100stmse  '\results_all_sbj.mat']);
                        case 2
                            load([path.auxres100 '\results_all_sbj.mat']);
                    end
            end
            
        % true positive only flag
        if TP_flag
            pval_CCA(find(DET_CCA ~= 1)) = NaN;
        end
        
        %% Calculate True/false positive/negative rates, precision, recall, ...
        tf_errors_AUX
    end
end

%% Average data
MSE = squeeze(nanmean(nanmean(nanmean(MSE_CCA,1),2),4));
CORR = squeeze(nanmean(nanmean(nanmean(CORR_CCA,1),2),4));
%FSCORE = squeeze(nanmean(nanmean(nanmean(CORR_CCA,1),2),4));



%% create imagesc matrices from data
% for ii=1:6
%     for jj=1:6
%         if 
%         else
%             CORRim(ii,jj)=NaN;
%         end
%     end
% end


mseffig = figure;
subplot(2,2,1)
imagesc

