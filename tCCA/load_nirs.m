load run1_0.1hz.nirs -mat

% auxilliary channels
acc1 = aux(:,2); % accelerometer
acc2 = aux(:,3); % accelerometer
acc3 = aux(:,4); % accelerometer
PPG = aux(:,5); % pulse
BP = aux(:,6); % blood pressure waveform
RESP = aux(:,7); % respiration

% get long distance channels (and active)
rhoSD_ssThresh = 15; % short distance threshold
ml = SD.MeasList;
mlAct = SD.MeasListAct;
lst = find(ml(:,4)==1);
rhoSD = zeros(length(lst),1);
posM = zeros(length(lst),3);
for iML = 1:length(lst)
    rhoSD(iML) = sum((SD.SrcPos(ml(lst(iML),1),:) - SD.DetPos(ml(lst(iML),2),:)).^2).^0.5;
    posM(iML,:) = (SD.SrcPos(ml(lst(iML),1),:) + SD.DetPos(ml(lst(iML),2),:)) / 2;
end
lstLong = lst(find(rhoSD>=rhoSD_ssThresh & mlAct(lst)==1));
lstShort = lst(find(rhoSD<=rhoSD_ssThresh & mlAct(lst)==1));

d_long = [d(:,lstLong), d(:,lstLong+size(d,2)/2)]; % first half 690 nm; second half 830 nm
d_short = [d(:,lstShort), d(:,lstShort+size(d,2)/2)]; % first half 690 nm; second half 830 nm