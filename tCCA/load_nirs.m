function [fq, t, AUX, d_long, d_short, d0_long, d0_short, d, d0, SD, s, lstLongAct,lstShortAct,lstHrfAdd] = load_nirs(filename)


if contains(filename,'nirs')
    load([filename], '-mat');
else
    load([filename '.nirs'], '-mat');
end
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


% get d_long and short(and active)
d_long = [d(:,lstLongAct), d(:,lstLongAct+size(d,2)/2)]; % first half 690 nm; second half 830 nm
d_short = [d(:,lstShortAct), d(:,lstShortAct+size(d,2)/2)]; % first half 690 nm; second half 830 nm
d0_long = [d0(:,lstLongAct), d0(:,lstLongAct+size(d,2)/2)]; % first half 690 nm; second half 830 nm
d0_short = [d0(:,lstShortAct), d0(:,lstShortAct+size(d,2)/2)]; % first half 690 nm; second half 830 nm