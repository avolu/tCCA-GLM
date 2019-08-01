%% Creates figure 4 (example average HRFs) of manuscript

% runCCA with subj 4, split 2, 
% parameter set: tlags = 3; stpsize = 2; cthresh = 0.3;
% regularized CCA on, plot_flag on.
%
% then run this script.

saveplot = true;

% plot channel
CH = 34;
% find index
idx=find(lstHrfAdd(:,1)==CH);

% time base
t = -2:1/fq:17;

figure
%% tCCA result
subplot(1,2,1)
plot(t(51:end),hrf.hrf_conc(1:426,1:2),'--k', 'LineWidth', 2)
hold on
plot(t,yavg_cca(:,1,CH),'r', 'LineWidth', 2)
plot(t,yavg_cca(:,2,CH),'b', 'LineWidth', 2)
grid on
xlim([-2,17])
ylim([-2.5,3.5]*1e-7)
xlabel('t / s')
ylabel('\Delta C / Mol')
title('50% HRF from GLM with CCA')

%% SS result
subplot(1,2,2)
plot(t(51:end),hrf.hrf_conc(1:426,1:2),'--k', 'LineWidth', 2)
hold on
plot(t,yavg_ss(:,1,CH),'r', 'LineWidth', 2)
plot(t,yavg_ss(:,2,CH),'b', 'LineWidth', 2)
grid on
xlim([-2,17])
ylim([-2.5,3.5]*1e-7)
xlabel('t / s')
ylabel('\Delta C / Mol')
title('50% HRF from GLM with SS')


%% display statistics
disp(['GLM with CCA - HbO/HbR. MSE: ' num2str(MSE_CCA(idx,1,2)) '/' num2str(MSE_CCA(idx,2,2)) ' Corr: ' num2str(CORR_CCA(idx,1,2)) '/' num2str(CORR_CCA(idx,2,2))])

disp(['GLM with SS - HbO/HbR. MSE: ' num2str(MSE_SS(idx,1,2)) '/' num2str(MSE_SS(idx,2,2)) ' Corr: ' num2str(CORR_SS(idx,1,2)) '/' num2str(CORR_SS(idx,2,2))])


set(gcf, 'Position',  [0,751,564,241])
path = 'C:\Users\avolu\Google Drive\tCCA_GLM_PAPER\FIGURES\Fig 4 hrf comp\';
if saveplot
    export_fig([path 'Fig4_hrf_comp'], '-pdf', '-transparent')
    export_fig([path 'Fig4_hrf_comp'], '-png', '-transparent', '-r300')
end