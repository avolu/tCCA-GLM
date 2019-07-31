clear all;

malexflag = 1;
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
filename = 'resting.nirs';
set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig))
set(groot,'defaultAxesCreateFcn',@(ax,~)set(ax.Toolbar,'Visible','off'))
sbjfolder = {'Subj33','Subj34','Subj36','Subj37','Subj38','Subj39', 'Subj40', 'Subj41', 'Subj43', 'Subj44','Subj46','Subj47','Subj49','Subj51'};


%% Options/Parameter Settings
rhoSD_ssThresh = 15;  % mm
flag_save = 0;
flag_conc = 1; % if 1 CCA inputs are in conc, if 0 CCA inputs are in intensity
segment = 100; % sec segment for fft
% results eval parameters
eval_param.HRFmin = -2;
eval_param.HRFmax = 17; % used only for block design runs
eval_param.Hb = 1; % 1 HbO / 0 HbR (for block only)
eval_param.pre = 5;  % HRF range in sec to calculate ttest
eval_param.post = 10;
% CCA parameters
flags.pcaf =  [0 0]; % no pca of X or AUX


% prune
flag_prune = false;
%plot flag
flag_plot = true;


% Optimum CCA parameters (obtained through cross-validation analysis of resting data -see tCCA folder)
timelag = 3; % in seconds
sts = 4 %  time step = 2, correct for 25 hz / 50 hz (delta t = 0.08s = 4/50) (this code is loading 50 hz data as opposed to rest of the functions)
cthresh = 0.3;


tic;


for sbj = 1:numel(sbjfolder) % loop across subjects
    disp(['subject #' num2str(sbj)]);
    
    
    %% load data
    load([path.dir filesep sbjfolder{sbj} filesep filename],'-mat');
    
    % sampling frequency
    fq = abs(1/(t(1)-t(2)));
    
    % auxilliary channels
    AUX = aux(:,2:7);
    % acc1 = aux(:,2); % accelerometer
    % acc2 = aux(:,3); % accelerometer
    % acc3 = aux(:,4); % accelerometer
    % PPG = aux(:,5); % pulse
    % BP = aux(:,6); % blood pressure waveform
    % RESP = aux(:,7); % respiration
    
    % lpf data entering CCA
    %     foo = d;
    foo = hmrBandpassFilt(d, fq, 0, 10);
    
    
    %% Prune noisy channels
    if flag_prune
        SD = enPruneChannels(d,SD,ones(size(d,1),1),[10000  10000000],5,[0  45],0);
    end
    
    %% Definition channel groups
    ml = SD.MeasList;
    lst = find(ml(:,4)==1); % 690
    rhoSD = zeros(length(lst),1);
    posM = zeros(length(lst),3);
    for iML = 1:length(lst)
        rhoSD(iML) = sum((SD.SrcPos(ml(lst(iML),1),:) - SD.DetPos(ml(lst(iML),2),:)).^2).^0.5;
        posM(iML,:) = (SD.SrcPos(ml(lst(iML),1),:) + SD.DetPos(ml(lst(iML),2),:)) / 2;
    end
    lstLongAct = lst(find(rhoSD> rhoSD_ssThresh & SD.MeasListAct(lst)==1)); % list of long and active channels
    lstShortAct = lst(find(rhoSD< rhoSD_ssThresh & SD.MeasListAct(lst)==1)); % list of long and active channels
    % get d_long and short(and active)
    d_long = [foo(:,lstLongAct), foo(:,lstLongAct+size(d,2)/2)]; % first half 690 nm; second half 830 nm
    d_short = [foo(:,lstShortAct), foo(:,lstShortAct+size(d,2)/2)]; % first half 690 nm; second half 830 nm
    
    
    %% lowpass filter AUX signals
    AUX = hmrBandpassFilt(AUX, fq, 0, 10);
    %% AUX signals
    AUX = [AUX, d_short]; % full AUX = [acc1 acc2 acc3 PPG BP RESP, d_short];
    %% zscore AUX signals
    AUX = zscore(AUX);
    
    mean_d_long = mean(d_long,1);
    
    
    %% Perform GLM with SS
    [yavg_ss, yavgstd_ss, tHRF, nTrialsSS, d_ss, yresid_ss, ysum2] = hmrDeconvTB_SS3rd(d, s, t, SD, [], [eval_param.HRFmin eval_param.HRFmax],0.5, 0.5, rhoSD_ssThresh);
    
    d_ss(:,[lstLongAct; lstLongAct+size(d,2)/2]) = d_ss(:,[lstLongAct; lstLongAct+size(d,2)/2])  +  repmat(mean_d_long, size(d_long,1),1);
    
    
    
    
    %% CCA EVAL
    param.tau = sts; %stepwidth for embedding in samples (tune to sample frequency!)
    param.NumOfEmb = ceil(timelag*fq / sts);
    param.ct = cthresh;   % correlation threshold
    
    X = d_long - repmat(mean_d_long, size(d_long,1),1);
    
    flags.shrink=true;
    [REG,  ADD] = rtcca(X,AUX,param,flags);
    
    
    %% Perform GLM with CCA
    [yavg_cca, yavgstd_cca, tHRF, nTrials, d_cca, yresid_cca, ysum2] = hmrDeconvTB_SS3rd(d, s, t, SD, REG, [eval_param.HRFmin eval_param.HRFmax],0.5, 0.5, 0);
    
    d_cca(:,[lstLongAct; lstLongAct+size(d,2)/2]) = d_cca(:,[lstLongAct; lstLongAct+size(d,2)/2])  +  repmat(mean_d_long, size(d_long,1),1);
    
    
    
    
    
    
    %% get power
    window = ones(1,fq*segment); %rectangular window
    
    for m = 1:size(lstLongAct);
        %% SS
        % 830
        [p_830,f]=pwelch(d_ss(:,lstLongAct(m) + size(d,2)/2),window,floor(segment*fq/2),0:1/segment:fq/2,fq);
        spectPow_830=sqrt(2*p_830/segment);  %duplicate due to both bands
        spectPow_830(1,:)=spectPow_830(1,:)/sqrt(2); %rectify DC
        spectPow_norm_830 = spectPow_830./spectPow_830(1); %normalize by DC
        % 690
        [p_690,f]=pwelch(d_ss(:,lstLongAct(m)),window,floor(segment*fq/2),0:1/segment:fq/2,fq);
        spectPow_690=sqrt(2*p_690/segment);  %duplicate due to both bands
        spectPow_690(1,:)=spectPow_690(1,:)/sqrt(2); %rectify DC
        spectPow_norm_690 = spectPow_690./spectPow_690(1); %normalize by DC
        
        % noise level
        f_lst_1 = find(f<=0.1 & f>=0.01);
        f_lst_2 = find(f<=0.5 & f>=0.1);
        f_lst_3 = find(f<=1.5 & f>=0.5);
        f_lst_4 = find(f<=10 & f>=1.5);
        f_lst_5 = find(f<=25 & f>=10);
        
        %830
        noise_d_ss_1_830(sbj,m) = mean(spectPow_norm_830(1,f_lst_1));
        noise_d_ss_2_830(sbj,m) = mean(spectPow_norm_830(1,f_lst_2));
        noise_d_ss_3_830(sbj,m) = mean(spectPow_norm_830(1,f_lst_3));
        noise_d_ss_4_830(sbj,m) = mean(spectPow_norm_830(1,f_lst_4));
        noise_d_ss_5_830(sbj,m) = mean(spectPow_norm_830(1,f_lst_5));
        
        %690
        noise_d_ss_1_690(sbj,m) = mean(spectPow_norm_690(1,f_lst_1));
        noise_d_ss_2_690(sbj,m) = mean(spectPow_norm_690(1,f_lst_2));
        noise_d_ss_3_690(sbj,m) = mean(spectPow_norm_690(1,f_lst_3));
        noise_d_ss_4_690(sbj,m) = mean(spectPow_norm_690(1,f_lst_4));
        noise_d_ss_5_690(sbj,m) = mean(spectPow_norm_690(1,f_lst_5));
        
        %% CCA
        % 830
        [p_830,f]=pwelch(d_cca(:,lstLongAct(m) + size(d,2)/2),window,floor(segment*fq/2),0:1/segment:fq/2,fq);
        spectPow_830=sqrt(2*p_830/segment);  %duplicate due to both bands
        spectPow_830(1,:)=spectPow_830(1,:)/sqrt(2); %rectify DC
        spectPow_norm_830 = spectPow_830./spectPow_830(1); %normalize by DC
        
        % 690
        [p_690,f]=pwelch(d_cca(:,lstLongAct(m)),window,floor(segment*fq/2),0:1/segment:fq/2,fq);
        spectPow_690=sqrt(2*p_690/segment);  %duplicate due to both bands
        spectPow_690(1,:)=spectPow_690(1,:)/sqrt(2); %rectify DC
        spectPow_norm_690 = spectPow_690./spectPow_690(1); %normalize by DC
        % noise level
        %830
        noise_d_cca_1_830(sbj,m) = mean(spectPow_norm_830(1,f_lst_1));
        noise_d_cca_2_830(sbj,m) = mean(spectPow_norm_830(1,f_lst_2));
        noise_d_cca_3_830(sbj,m) = mean(spectPow_norm_830(1,f_lst_3));
        noise_d_cca_4_830(sbj,m) = mean(spectPow_norm_830(1,f_lst_4));
        noise_d_cca_5_830(sbj,m) = mean(spectPow_norm_830(1,f_lst_5));
        
        %690
        noise_d_cca_1_690(sbj,m) = mean(spectPow_norm_690(1,f_lst_1));
        noise_d_cca_2_690(sbj,m) = mean(spectPow_norm_690(1,f_lst_2));
        noise_d_cca_3_690(sbj,m) = mean(spectPow_norm_690(1,f_lst_3));
        noise_d_cca_4_690(sbj,m) = mean(spectPow_norm_690(1,f_lst_4));
        noise_d_cca_5_690(sbj,m) = mean(spectPow_norm_690(1,f_lst_5));
        
        
        
    end
    
    
    if flag_plot
        figure
        for m = 1:size(lstLongAct);
            %% SS
            [p,f]=pwelch(d_ss(:,lstLongAct(m) + size(d,2)/2),window,floor(segment*fq/2),0:1/segment:fq/2,fq);
            spectPow=sqrt(2*p/segment);  %duplicate due to both bands
            spectPow(1,:)=spectPow(1,:)/sqrt(2); %rectify DC
            spectPow_norm = spectPow./spectPow(1); %normalize by DC
            spectPow_norm_ss(:,m) = spectPow_norm;
            h1 = loglog(f,spectPow_norm,'b','LineWidth',0.1);
            h1.Color=[0,0,1,0.4];
            hold on;
            
            %% CCA
            [p,f]=pwelch(d_cca(:,lstLongAct(m) + size(d,2)/2),window,floor(segment*fq/2),0:1/segment:fq/2,fq);
            spectPow=sqrt(2*p/segment);  %duplicate due to both bands
            spectPow(1,:)=spectPow(1,:)/sqrt(2); %rectify DC
            spectPow_norm = spectPow./spectPow(1); %normalize by DC
            spectPow_norm_cca(:,m) = spectPow_norm;
            h2 = loglog(f,spectPow_norm,'r','LineWidth',0.1);
            h2.Color=[0.9,0,0,0.4];
            
            
        end
        legend('\color{blue} GLM with SS', '\color{red} GLM with CCA');
        
        
        % plot average of all active channels
        %%  ss
        p_avg = mean(spectPow_norm_ss,2); %    p_avg = mean(spectPow_norm_ss(:,lstLongAct),2);
        loglog(f,p_avg,'b','LineWidth',2);
        
        %% cca
        p_avg = mean(spectPow_norm_cca,2);
        loglog(f,p_avg,'r','LineWidth',2);
        grid on;
        xlim([min(f) max(f)]);
        ylim([1e-6 0.1]);
        xlabel('Frequency [Hz]');
        ylabel('d830 power');
        
        
        
        
        figure;
        p_avg = mean(spectPow_norm_ss,2); %    p_avg = mean(spectPow_norm_ss(:,lstLongAct),2);
        loglog(f,p_avg,'b','LineWidth',2);
        hold on;
        %% cca
        p_avg = mean(spectPow_norm_cca,2);
        loglog(f,p_avg,'r','LineWidth',2);
        grid on;
        xlim([min(f) max(f)]);
        %                     ylim([min(min([spectPow_norm_ss(:) spectPow_norm_cca(:)])) max(max([spectPow_norm_ss(:) spectPow_norm_cca(:)]))]);
        ylim([1e-5 0.1]);
        
        xlabel('Frequency [Hz]');
        ylabel('d830 power');
        legend('\color{blue} GLM with SS', '\color{red} GLM with CCA');
    end
    
    
end

%% Figure 12 in paper
%% plot power at different frequency bands *for 830*

figure;
subplot(2,2,1);
scatter(reshape(noise_d_ss_1_830,1,size(noise_d_ss_1_830,1)*size(noise_d_ss_1_830,2)),reshape(noise_d_cca_1_830,1,size(noise_d_cca_1_830,1)*size(noise_d_cca_1_830,2)),'bo','filled','MarkerFaceAlpha',3/8);hold on;
plot(mean(noise_d_ss_1_830(:)),mean(noise_d_cca_1_830(:)),'bp','MarkerSize',15,'MarkerEdgeColor','k','MarkerFaceColor','b','LineWidth',2);
plot([1e-6 0.1], [1e-6 0.1],'k');
set(gca,'xscale','log')
set(gca,'yscale','log')
legend('\color{blue} 0.01 < f < 0.1 Hz');
xlabel('Power- GLM with SS');
ylabel('Power - GLM with CCA');
xlim([min(min([noise_d_ss_1_830(:) noise_d_cca_1_830(:)])),max(max([noise_d_ss_1_830(:) noise_d_cca_1_830(:)]))])
ylim([min(min([noise_d_ss_1_830(:) noise_d_cca_1_830(:)])),max(max([noise_d_ss_1_830(:) noise_d_cca_1_830(:)]))])
grid;

subplot(2,2,2);
scatter(reshape(noise_d_ss_2_830,1,size(noise_d_ss_1_830,1)*size(noise_d_ss_1_830,2)),reshape(noise_d_cca_2_830,1,size(noise_d_cca_1_830,1)*size(noise_d_cca_1_830,2)),'go','filled','MarkerFaceAlpha',3/8);hold on;
plot(mean(noise_d_ss_2_830(:)),mean(noise_d_cca_2_830(:)),'gp','MarkerSize',15,'MarkerEdgeColor','k','MarkerFaceColor','g','LineWidth',2);
plot([1e-6 0.1], [1e-6 0.1],'k');
set(gca,'xscale','log')
set(gca,'yscale','log')
legend('\color{green} 0.1 < f < 0.5 Hz');
xlabel('Power- GLM with SS');
ylabel('Power - GLM with CCA');
xlim([min(min([noise_d_ss_2_830(:) noise_d_cca_2_830(:)])),max(max([noise_d_ss_2_830(:) noise_d_cca_2_830(:)]))])
ylim([min(min([noise_d_ss_2_830(:) noise_d_cca_2_830(:)])),max(max([noise_d_ss_2_830(:) noise_d_cca_2_830(:)]))])
grid;

subplot(2,2,3);
scatter(reshape(noise_d_ss_3_830,1,size(noise_d_ss_1_830,1)*size(noise_d_ss_1_830,2)),reshape(noise_d_cca_3_830,1,size(noise_d_cca_1_830,1)*size(noise_d_cca_1_830,2)),'ro','filled','MarkerFaceAlpha',3/8);hold on;
plot(mean(noise_d_ss_3_830(:)),mean(noise_d_cca_3_830(:)),'rp','MarkerSize',15,'MarkerEdgeColor','k','MarkerFaceColor','r','LineWidth',2);
plot([1e-6 0.1], [1e-6 0.1],'k');
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('Power- GLM with SS');
ylabel('Power - GLM with CCA');
xlim([min(min([noise_d_ss_3_830(:) noise_d_cca_1_830(:)])),max(max([noise_d_ss_3_830(:) noise_d_cca_1_830(:)]))])
ylim([min(min([noise_d_ss_3_830(:) noise_d_cca_1_830(:)])),max(max([noise_d_ss_3_830(:) noise_d_cca_1_830(:)]))])
legend('\color{red} 0.5 < f < 1.5 Hz')
grid;

subplot(2,2,4);
scatter(reshape(noise_d_ss_4_830,1,size(noise_d_ss_1_830,1)*size(noise_d_ss_1_830,2)),reshape(noise_d_cca_4_830,1,size(noise_d_cca_1_830,1)*size(noise_d_cca_1_830,2)),'mo','filled','MarkerFaceAlpha',3/8);hold on;
plot(mean(noise_d_ss_4_830(:)),mean(noise_d_cca_4_830(:)),'mp','MarkerSize',15,'MarkerEdgeColor','k','MarkerFaceColor','m','LineWidth',2);
plot([1e-6 0.1], [1e-6 0.1],'k');
legend('\color{magenta} 1.5 < f < 10 Hz');
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('Power- GLM with SS');
ylabel('Power - GLM with CCA');
xlim([min(min([noise_d_ss_4_830(:) noise_d_cca_4_830(:)])),max(max([noise_d_ss_4_830(:) noise_d_cca_4_830(:)]))])
ylim([min(min([noise_d_ss_4_830(:) noise_d_cca_4_830(:)])),max(max([noise_d_ss_4_830(:) noise_d_cca_4_830(:)]))])
grid;


% plot bar plot for mean and std
noise_ss= [mean(noise_d_ss_1_830(:)),mean(noise_d_ss_2_830(:)),mean(noise_d_ss_3_830(:)),mean(noise_d_ss_4_830(:))]; % SS
noise_cca= [mean(noise_d_cca_1_830(:)),mean(noise_d_cca_2_830(:)),mean(noise_d_cca_3_830(:)),mean(noise_d_cca_4_830(:))]; % SS
[h,p_noise_1,c,stats]=ttest(noise_d_ss_1_830(:),noise_d_cca_1_830(:));
fprintf('p_value for t-test between SS and tCCA for 0.01 < f < 0.1 Hz:   %.4d \n',p_noise_1);
[h,p_noise_2,c,stats]=ttest(noise_d_ss_2_830(:),noise_d_cca_2_830(:));
fprintf('p_value for t-test between SS and tCCA for 0.1 < f < 0.5 Hz:   %.4d \n',p_noise_2);
[h,p_noise_3,c,stats]=ttest(noise_d_ss_3_830(:),noise_d_cca_3_830(:));
fprintf('p_value for t-test between SS and tCCA for 0.5 < f < 1.5 Hz:   %.4d \n',p_noise_3);
[h,p_noise_4,c,stats]=ttest(noise_d_ss_4_830(:),noise_d_cca_4_830(:));
fprintf('p_value for t-test between SS and tCCA for 1.5 < f < 10 Hz:   %.4d \n',p_noise_4);

group = [noise_ss',noise_cca'];
foo= [noise_ss',noise_cca']'; foo = foo(:)';

SEM=[std(noise_d_ss_1_830(:)),std(noise_d_cca_1_830(:)),std(noise_d_ss_2_830(:)),std(noise_d_cca_2_830(:)),std(noise_d_ss_3_830(:)),std(noise_d_cca_3_830(:)),std(noise_d_ss_4_830(:)),std(noise_d_cca_4_830(:))];  % values for error bars
SEM = SEM./sqrt(size(noise_d_ss_1_830(:),1)-1);

figure
hold on
bar(1:4,group)
x_shift = 0.14;
x = [1-x_shift,1+x_shift,2-x_shift,2+x_shift,3-x_shift,3+x_shift,4-x_shift,4+x_shift];
errorbar(x,foo,SEM,'.') %errorbar(x,y,err)
xticklabels({'0.01 < f < 0.1 Hz',' ','0.1 < f < 0.5 Hz',' ','0.5 < f < 1.5 Hz',' ','1.5 < f < 10 Hz'});
xtickangle(30);
legend('\color{blue} GLM with SS', '\color{red} GLM with CCA');
title('830')



%% plot power at different frequency bands *for 690*

figure;
subplot(2,2,1);
scatter(reshape(noise_d_ss_1_690,1,size(noise_d_ss_1_690,1)*size(noise_d_ss_1_690,2)),reshape(noise_d_cca_1_690,1,size(noise_d_cca_1_690,1)*size(noise_d_cca_1_690,2)),'bo','filled','MarkerFaceAlpha',3/8);hold on;
plot(mean(noise_d_ss_1_690(:)),mean(noise_d_cca_1_690(:)),'bp','MarkerSize',15,'MarkerEdgeColor','k','MarkerFaceColor','b','LineWidth',2);
plot([1e-6 0.1], [1e-6 0.1],'k');
set(gca,'xscale','log')
set(gca,'yscale','log')
legend('\color{blue} 0.01 < f < 0.1 Hz');
xlabel('Power- GLM with SS');
ylabel('Power - GLM with CCA');
xlim([min(min([noise_d_ss_1_690(:) noise_d_cca_1_690(:)])),max(max([noise_d_ss_1_690(:) noise_d_cca_1_690(:)]))])
ylim([min(min([noise_d_ss_1_690(:) noise_d_cca_1_690(:)])),max(max([noise_d_ss_1_690(:) noise_d_cca_1_690(:)]))])
grid;

subplot(2,2,2);
scatter(reshape(noise_d_ss_2_690,1,size(noise_d_ss_1_690,1)*size(noise_d_ss_1_690,2)),reshape(noise_d_cca_2_690,1,size(noise_d_cca_1_690,1)*size(noise_d_cca_1_690,2)),'go','filled','MarkerFaceAlpha',3/8);hold on;
plot(mean(noise_d_ss_2_690(:)),mean(noise_d_cca_2_690(:)),'gp','MarkerSize',15,'MarkerEdgeColor','k','MarkerFaceColor','g','LineWidth',2);
plot([1e-6 0.1], [1e-6 0.1],'k');
set(gca,'xscale','log')
set(gca,'yscale','log')
legend('\color{green} 0.1 < f < 0.5 Hz');
xlabel('Power- GLM with SS');
ylabel('Power - GLM with CCA');
xlim([min(min([noise_d_ss_2_690(:) noise_d_cca_2_690(:)])),max(max([noise_d_ss_2_690(:) noise_d_cca_2_690(:)]))])
ylim([min(min([noise_d_ss_2_690(:) noise_d_cca_2_690(:)])),max(max([noise_d_ss_2_690(:) noise_d_cca_2_690(:)]))])
grid;

subplot(2,2,3);
scatter(reshape(noise_d_ss_3_690,1,size(noise_d_ss_1_690,1)*size(noise_d_ss_1_690,2)),reshape(noise_d_cca_3_690,1,size(noise_d_cca_1_690,1)*size(noise_d_cca_1_690,2)),'ro','filled','MarkerFaceAlpha',3/8);hold on;
plot(mean(noise_d_ss_3_690(:)),mean(noise_d_cca_3_690(:)),'rp','MarkerSize',15,'MarkerEdgeColor','k','MarkerFaceColor','r','LineWidth',2);
plot([1e-6 0.1], [1e-6 0.1],'k');
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('Power- GLM with SS');
ylabel('Power - GLM with CCA');
xlim([min(min([noise_d_ss_3_690(:) noise_d_cca_1_690(:)])),max(max([noise_d_ss_3_690(:) noise_d_cca_1_690(:)]))])
ylim([min(min([noise_d_ss_3_690(:) noise_d_cca_1_690(:)])),max(max([noise_d_ss_3_690(:) noise_d_cca_1_690(:)]))])
legend('\color{red} 0.5 < f < 1.5 Hz')
grid;

subplot(2,2,4);
scatter(reshape(noise_d_ss_4_690,1,size(noise_d_ss_1_690,1)*size(noise_d_ss_1_690,2)),reshape(noise_d_cca_4_690,1,size(noise_d_cca_1_690,1)*size(noise_d_cca_1_690,2)),'mo','filled','MarkerFaceAlpha',3/8);hold on;
plot(mean(noise_d_ss_4_690(:)),mean(noise_d_cca_4_690(:)),'mp','MarkerSize',15,'MarkerEdgeColor','k','MarkerFaceColor','m','LineWidth',2);
plot([1e-6 0.1], [1e-6 0.1],'k');
legend('\color{magenta} 1.5 < f < 10 Hz');
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('Power- GLM with SS');
ylabel('Power - GLM with CCA');
xlim([min(min([noise_d_ss_4_690(:) noise_d_cca_4_690(:)])),max(max([noise_d_ss_4_690(:) noise_d_cca_4_690(:)]))])
ylim([min(min([noise_d_ss_4_690(:) noise_d_cca_4_690(:)])),max(max([noise_d_ss_4_690(:) noise_d_cca_4_690(:)]))])
grid;


% plot bar plot for mean and std
noise_ss= [mean(noise_d_ss_1_690(:)),mean(noise_d_ss_2_690(:)),mean(noise_d_ss_3_690(:)),mean(noise_d_ss_4_690(:))]; % SS
noise_cca= [mean(noise_d_cca_1_690(:)),mean(noise_d_cca_2_690(:)),mean(noise_d_cca_3_690(:)),mean(noise_d_cca_4_690(:))]; % SS
[h,p_noise_1,c,stats]=ttest(noise_d_ss_1_690(:),noise_d_cca_1_690(:));
[h,p_noise_2,c,stats]=ttest(noise_d_ss_2_690(:),noise_d_cca_2_690(:));
[h,p_noise_3,c,stats]=ttest(noise_d_ss_3_690(:),noise_d_cca_3_690(:));
[h,p_noise_4,c,stats]=ttest(noise_d_ss_4_690(:),noise_d_cca_4_690(:));

group = [noise_ss',noise_cca'];
foo= [noise_ss',noise_cca']'; foo = foo(:)';

SEM=[std(noise_d_ss_1_690(:)),std(noise_d_cca_1_690(:)),std(noise_d_ss_2_690(:)),std(noise_d_cca_2_690(:)),std(noise_d_ss_3_690(:)),std(noise_d_cca_3_690(:)),std(noise_d_ss_4_690(:)),std(noise_d_cca_4_690(:))];  % values for error bars
SEM = SEM./sqrt(size(noise_d_ss_1_690(:),1)-1);

figure
hold on
bar(1:4,group)
x_shift = 0.14;
x = [1-x_shift,1+x_shift,2-x_shift,2+x_shift,3-x_shift,3+x_shift,4-x_shift,4+x_shift];
errorbar(x,foo,SEM,'.') %errorbar(x,y,err)
xticklabels({'0.01 < f < 0.1 Hz',' ','0.1 < f < 0.5 Hz',' ','0.5 < f < 1.5 Hz',' ','1.5 < f < 10 Hz'});
xtickangle(30);
legend('\color{blue} GLM with SS', '\color{red} GLM with CCA');
title('690')
