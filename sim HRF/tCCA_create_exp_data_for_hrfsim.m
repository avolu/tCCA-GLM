path.dir = 'C:\Users\mayucel\Google Drive\CCA_PAPER\';
path.fname = 'resting.nirs';
path.hrfname = 'hrf_simdat.mat';
path.savename = 'resting_sim.nirs';
addpath(genpath([path.dir filesep 'CODE']));
sbjfolder = {'Subj33','Subj34','Subj36','Subj37','Subj38','Subj39', 'Subj40', 'Subj41', 'Subj43', 'Subj44','Subj46','Subj47','Subj49','Subj51'};


for ss = 1:numel(sbjfolder) % loop across subjects
% load participant data
nirs = load([path.dir 'FB_RESTING_DATA' filesep sbjfolder{ss} filesep path.fname], '-mat');
% load simulated hrf
hrf = load([path.dir 'CODE' filesep 'sim HRF' filesep path.hrfname]);
% hrf = load([path.dir path.hrfname]);

%% call the function
[nirs_hrf] = addSimHRF(nirs, hrf, false);

save([path.dir 'FB_RESTING_DATA' filesep sbjfolder{ss} filesep path.savename], '-struct', 'nirs_hrf')
disp('data saved.')
end





