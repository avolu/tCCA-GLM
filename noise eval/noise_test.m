
%% Create simulated data and parameters
fs = 1 / 0.01;
t=1:1/fs:100;


%% Sources
% LFO
Flfo=0.1;
Slfo = sin(2*pi*Flfo*t)';
% Heart rate
Fhr=1;
Shr = sin(2*pi*Fhr*t)';
% Breathing
Fbr=0.3;
Sbr = sin(2*pi*Fbr*t)';
% white noise with -3dBW / 0.5W
Snoise = wgn(numel(t),1,-3);
% assemble:
S = [Slfo Shr Sbr Snoise];
%% Mixing matrix
A = rand(4,4);
%% Measurement channels X
X = S*A;

%% "Regressor"
Pn = 0;
Nfac = 1;
rnoise = wgn(numel(t),1,Pn);
Rhr = Shr + Nfac*rnoise;

%% Remove noisy heart signal from channels X with perfect known scaling factor
X_clean = X - [Rhr*A(2,1) Rhr*A(2,2) Rhr*A(2,3) Rhr*A(2,4)];

%% Plot 
xl = [0 25];
% S
tlab = {'Source LFO', 'Source Heart', 'Source Breating', 'GWN, -3dbW'}
figure
for ii=1:4
subplot(2,2,ii)
plot(t,S(:,ii))
xlim(xl)
xlabel('t/s')
title(tlab{ii})
end

% X
figure
for ii=1:4
subplot(2,2,ii)
plot(t,X(:,ii))
xlim(xl)
xlabel('t/s')
title(['mixed signals CH' num2str(ii)])
end

% Mixed signals and corresponding weighted Regressor
figure
for ii=1:4
subplot(2,2,ii)
plot(t,X(:,ii))
hold on
plot(t,Rhr*A(2,ii))
xlim(xl)
xlabel('t/s')
title(['mixed CH' num2str(ii) ' & regressor, weighted with ' num2str(A(2,ii))])
end

% Cleaned X
figure
for ii=1:4
subplot(2,2,ii)
plot(t,X_clean(:,ii))
xlim(xl)
xlabel('t/s')
title(['cleaned signals CH' num2str(ii)])
end

% Spectra before and after
figure
for ii=1:4
subplot(2,2,ii)
[px,f] = pspectrum(X(:,ii),fs);
[pxc,f] = pspectrum(X_clean(:,ii),fs);
loglog(f,px)
hold on
loglog(f,pxc)
xlabel('f/Hz')
axis tight
title(['cleaned signals CH' num2str(ii)])
end
