% load visual_probe_plot.SD -mat;
%% plot cross-trial Mean with STD
function plot_block(MEAN_null, MEAN_new, HRFmin, HRFmax, fq, pOxy_null, pOxy_new,ss,STD_null,STD_new, tHRF, timelag,fpath,lstHrfAdd)
cd(fpath)
load visual_probe_plot.SD -mat
coor_xy = [(SD.SrcPos(SD.MeasList(1:size(SD.MeasList,1)/2,1),1) + SD.DetPos(SD.MeasList(1:size(SD.MeasList,1)/2,2),1))/5,(SD.SrcPos(SD.MeasList(1:size(SD.MeasList,1)/2,1),2)+ SD.DetPos(SD.MeasList(1:size(SD.MeasList,1)/2,2),2))/5];
% F=size(K);
a = coor_xy(:,1);
b = coor_xy(:,2);

nf = 90;
% baseline correct the mean for the plots
Baseline = mean(MEAN_null(1:abs(fq*HRFmin),:),1);
MEAN_null = MEAN_null - Baseline;
Baseline = mean(MEAN_new(1:abs(fq*HRFmin),:),1);
MEAN_new = MEAN_new - Baseline;


MEAN_null = downsample(MEAN_null,nf);
STD_null = downsample(STD_null,nf);
MEAN_new = downsample(MEAN_new,nf);
STD_new = downsample(STD_new,nf);
tHRF = downsample(tHRF,nf);


% % z1 = 0.02; % normalize
% % z2 = 0.95;

z1 = 0.12; % normalize
z2 = 0.75;% adjust this for different probes

a1 = (a*(z2-z1)/(max(a)-min(a)));
a1 = a1+ z1-min(a1);

b1 = ( b*(z2-z1)/(max(b)-min(b)));
b1 = b1+ z1-min(b1);

%get rid of the short-separation channels
rhoSD_ssThresh = 15;
ml = SD.MeasList;
mlAct = SD.MeasListAct;
lst = find(ml(:,4)==1);
rhoSD = zeros(length(lst),1);
posM = zeros(length(lst),3);
for iML = 1:length(lst)
    rhoSD(iML) = sum((SD.SrcPos(ml(lst(iML),1),:) - SD.DetPos(ml(lst(iML),2),:)).^2).^0.5;
    posM(iML,:) = (SD.SrcPos(ml(lst(iML),1),:) + SD.DetPos(ml(lst(iML),2),:)) / 2;
end
lstLL = lst(find(rhoSD>=rhoSD_ssThresh & mlAct(lst)==1));

ylim1 = -1e-6;
ylim2 = 1e-6;
xlim1 = HRFmin;
xlim2 = HRFmax;
figure;
for i =lstLL'
    h=subplot('Position',[a1(i),b1(i),0.06,0.1]);
    if pOxy_null(ss,i)<=0.05
        errorbar(min(tHRF):(max(tHRF) -min(tHRF))/(size(MEAN_null,1)-1):max(tHRF),MEAN_null(:,i),STD_null(:,i),STD_null(:,i),'r','LineWidth',2);
        title(['    p = ' (num2str(pOxy_null(ss,i),1))],'FontSize',15,'FontWeight','bold','color','k') ;
    elseif pOxy_null(ss,i)>0.05
        errorbar(min(tHRF):(max(tHRF) -min(tHRF))/(size(MEAN_null,1)-1):max(tHRF),MEAN_null(:,i),STD_null(:,i),STD_null(:,i),'color',[0.5 0.5 0.5],'LineWidth',2);
    end
    if any(lstHrfAdd(:,1) == i)
        xlabel('HRF')
    end
    grid; ylim([ylim1 ylim2]);xlim([xlim1 xlim2]);
end
suptitle(['GLM - Subject # ' num2str(ss) ' sec']) ;


figure;
for i=lstLL'
    h=subplot('Position',[a1(i),b1(i),0.06,0.1]);
    if pOxy_new(ss,i)<=0.05
        errorbar(min(tHRF):(max(tHRF) -min(tHRF))/(size(MEAN_new,1)-1):max(tHRF),MEAN_new(:,i),STD_new(:,i),STD_new(:,i),'r','LineWidth',2);
        title(['    p = ' (num2str(pOxy_new(ss,i),1))],'FontSize',15,'FontWeight','bold','color','k') ;
    elseif pOxy_new(ss,i)>0.05
        errorbar(min(tHRF):(max(tHRF) -min(tHRF))/(size(MEAN_new,1)-1):max(tHRF),MEAN_new(:,i),STD_new(:,i),STD_new(:,i),'color',[0.5 0.5 0.5],'LineWidth',2);
    end
    if any(lstHrfAdd(:,1) == i)
        xlabel('HRF')
    end
    grid; ylim([ylim1 ylim2]); xlim([xlim1 xlim2]);
end
suptitle(['GLM with CCA - Subject # ' num2str(ss) ',  t_l_a_g= ' num2str(timelag) ' sec']) ;
