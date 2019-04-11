%% EXAMPLE SCRIPT
% BE AWARE:
% filtering and zscoring might be necessary/helpful. CCA components
% (regressors) have arbitrary scaling. Should not be a problem in the GLM
% regression, but keep it in mind!

% load data
load_nirs

% merge data to right format
AUX = [acc1 acc2 acc3 BP RESP aux d_short];
X = d_long;
fs=50; %sample frequency

% Lowpass-Filter data?
filtflag = true

if filtflag
    fclow = 20;
    order = 3;
    [f2,f1] = butter(order, fclow/fs*2, 'low');
    X = filtfilt(f2,f1,X);
    AUX = filtfilt(f2,f1,AUX);
end


param.tau= 2% fsmpl/5;   %stepwidth for embedding in samples (tune to sample frequency!)
param.NumOfEmb= 75; %number of total embeddings (total time window will be tau * NumOfEmb)
param.ct= 0.5;   % correlation threshold
flags.pcaf=  [0 0]; % no pca of X or AUX

%% Perform the shiny script
[REG, ADD] = perf_temp_emb_cca(X,AUX,param,flags);

%% Plot regressors
figure
plot(REG)

%% FFT
figure;
pspectrum(REG, fs);
% alternatively:
figure;
segment = 100; % sec
foo = REG(:,1);
[p,f]=pwelch(foo,ones(1,fs*segment),floor(segment*fs/2),0:1/segment:fs/2,fs);
spectPow=sqrt(2*p/segment);  %duplicate due to both bands
spectPow(1,:)=spectPow(1,:)/sqrt(2); %rectify DC
semilogy(f,spectPow);
xlim([0 25]);
xlabel('Frequency [Hz]');
ylabel('XXX power');


