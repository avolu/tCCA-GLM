% load visual_probe_plot.SD -mat;
%% plot cross-trial Mean with STD
function plot_block(MEAN_SS, MEAN_CCA, CORR_SS, CORR_CCA, MSE_SS, MSE_CCA, HRFmin, HRFmax, fq, pOxy_SS, pOxy_CCA,ss,STD_SS,STD_CCA, tHRF, timelag,fpath,lstHrfAdd,hrf)

cd(fpath)
load visual_probe_plot.SD -mat
coor_xy = [(SD.SrcPos(SD.MeasList(1:size(SD.MeasList,1)/2,1),1) + SD.DetPos(SD.MeasList(1:size(SD.MeasList,1)/2,2),1))/5,(SD.SrcPos(SD.MeasList(1:size(SD.MeasList,1)/2,1),2)+ SD.DetPos(SD.MeasList(1:size(SD.MeasList,1)/2,2),2))/5];
% F=size(K);
a = coor_xy(:,1);
b = coor_xy(:,2);

nf = 40;


% baseline correct the mean for the plots
% for HbO and HbR
for i = 1:size(MEAN_SS,3)
    for j = 1:size(MEAN_SS,2)
        Baseline(j,i) = squeeze(mean(MEAN_SS(1:abs(fq*HRFmin),j,i),1));
        MEAN_SS(:,j,i) = squeeze(MEAN_SS(:,j,i)) - squeeze(Baseline(j,i));
        Baseline(j,i) = mean(MEAN_CCA(1:abs(fq*HRFmin),j,i),1);
        MEAN_CCA(:,j,i) = MEAN_CCA(:,j,i) - Baseline(j,i);
    end
end

% downsample for visualization
for i = 1:size(MEAN_SS,3)
    for j = 1:size(MEAN_SS,2)
        MEAN_SS_down(:,j,i) = downsample(MEAN_SS(:,j,i),nf);
        STD_SS_down(:,j,i) = downsample(STD_SS(:,j,i),nf);
        MEAN_CCA_down(:,j,i) = downsample(MEAN_CCA(:,j,i),nf);
        STD_CCA_down(:,j,i) = downsample(STD_CCA(:,j,i),nf);
    end
end
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
ylim2 = 1.5e-6;
xlim1 = HRFmin;
xlim2 = HRFmax;


figure;
% SS
j = 1; % HbO
foo = 1;
for i =lstLL'
    h=subplot('Position',[a1(i),b1(i),0.06,0.1]);
    hold on;
    if pOxy_SS(i,j)<=0.05
        errorbar(min(tHRF):(max(tHRF) -min(tHRF))/(size(MEAN_SS_down,1)-1):max(tHRF),MEAN_SS_down(:,i,j),STD_SS_down(:,i,j),STD_SS_down(:,i,j),'r','LineWidth',2);
        title(['    p = ' (num2str(pOxy_SS(i,j),1))],'FontSize',15,'FontWeight','bold','color','k') ;
    elseif pOxy_SS(i,j)>0.05
        errorbar(min(tHRF):(max(tHRF) -min(tHRF))/(size(MEAN_SS_down,1)-1):max(tHRF),MEAN_SS_down(:,i,j),STD_SS_down(:,i,j),STD_SS_down(:,i,j),'color',[0.5 0.5 0.5],'LineWidth',2);
    end
    if any(lstHrfAdd(:,1) == i)
        First_line = ['HRF, Corr: ' num2str(CORR_SS(foo,j),'%0.2g')];
        Second_line =  [' MSE: ' num2str(MSE_SS(foo,j),'%0.2g') ];
        xlabel({First_line;Second_line})
        foo = foo + 1;
        plot([0:1/fq:max(hrf.t_hrf+1/fq)],hrf.hrf_conc(:,j));
    end
    txt = ['ch ' num2str(i)];
    ylabel(txt);
    grid; ylim([ylim1 ylim2]);xlim([xlim1 xlim2]);
end
suptitle(['GLM with SS - Subject # ' num2str(ss) ' sec']) ;


figure;
% CCA
foo = 1;
for i =lstLL'
    j = 1; % HbO
    h=subplot('Position',[a1(i),b1(i),0.06,0.1]);
    hold on;
    if pOxy_CCA(i,j)<=0.05
        errorbar(min(tHRF):(max(tHRF) -min(tHRF))/(size(MEAN_CCA_down,1)-1):max(tHRF),MEAN_CCA_down(:,i,j),STD_CCA_down(:,i,j),STD_CCA_down(:,i,j),'r','LineWidth',2);
        title(['    p = ' (num2str(pOxy_CCA(i,j),1))],'FontSize',15,'FontWeight','bold','color','k') ;
    elseif pOxy_CCA(i,j)>0.05
        errorbar(min(tHRF):(max(tHRF) -min(tHRF))/(size(MEAN_CCA_down,1)-1):max(tHRF),MEAN_CCA_down(:,i,j),STD_CCA_down(:,i,j),STD_CCA_down(:,i,j),'color',[0.5 0.5 0.5],'LineWidth',2);
    end
    if any(lstHrfAdd(:,1) == i)
        First_line = ['HRF, Corr: ' num2str(CORR_CCA(foo,j),'%0.2g')];
        Second_line =  [' MSE: ' num2str(MSE_CCA(foo,j),'%0.2g') ];
        xlabel({First_line;Second_line})
        foo = foo + 1;
        plot([0:1/fq:max(hrf.t_hrf+1/fq)],hrf.hrf_conc(:,j));
    end
    txt = ['ch ' num2str(i)];
    ylabel(txt);
    grid; ylim([ylim1 ylim2]);xlim([xlim1 xlim2]);
end
suptitle(['GLM with CCA - Subject # ' num2str(ss) ',  t_l_a_g= ' num2str(timelag) ' sec']) ;







figure;
% SS
j = 2; % HbR
foo = 1;
for i =lstLL'
    h=subplot('Position',[a1(i),b1(i),0.06,0.1]);
    hold on;
    if pOxy_SS(i,j)<=0.05
        errorbar(min(tHRF):(max(tHRF) -min(tHRF))/(size(MEAN_SS_down,1)-1):max(tHRF),MEAN_SS_down(:,i,j),STD_SS_down(:,i,j),STD_SS_down(:,i,j),'r','LineWidth',2);
        title(['    p = ' (num2str(pOxy_SS(i,j),1))],'FontSize',15,'FontWeight','bold','color','k') ;
    elseif pOxy_SS(i,j)>0.05
        errorbar(min(tHRF):(max(tHRF) -min(tHRF))/(size(MEAN_SS_down,1)-1):max(tHRF),MEAN_SS_down(:,i,j),STD_SS_down(:,i,j),STD_SS_down(:,i,j),'color',[0.5 0.5 0.5],'LineWidth',2);
    end
    if any(lstHrfAdd(:,1) == i)
        First_line = ['HRF, Corr: ' num2str(CORR_SS(foo,j),'%0.2g')];
        Second_line =  [' MSE: ' num2str(MSE_SS(foo,j),'%0.2g') ];
        xlabel({First_line;Second_line})
        foo = foo + 1;
        plot([0:1/fq:max(hrf.t_hrf+1/fq)],hrf.hrf_conc(:,j));
    end
    txt = ['ch ' num2str(i)];
    ylabel(txt);
    grid; ylim([ylim1 ylim2]);xlim([xlim1 xlim2]);
end
suptitle(['GLM with SS - Subject # ' num2str(ss) ' sec']) ;


figure;
% CCA
foo = 1;
for i =lstLL'
    j = 2; % HbR
    h=subplot('Position',[a1(i),b1(i),0.06,0.1]);
    hold on;
    if pOxy_CCA(i,j)<=0.05
        errorbar(min(tHRF):(max(tHRF) -min(tHRF))/(size(MEAN_CCA_down,1)-1):max(tHRF),MEAN_CCA_down(:,i,j),STD_CCA_down(:,i,j),STD_CCA_down(:,i,j),'r','LineWidth',2);
        title(['    p = ' (num2str(pOxy_CCA(i,j),1))],'FontSize',15,'FontWeight','bold','color','k') ;
    elseif pOxy_CCA(i,j)>0.05
        errorbar(min(tHRF):(max(tHRF) -min(tHRF))/(size(MEAN_CCA_down,1)-1):max(tHRF),MEAN_CCA_down(:,i,j),STD_CCA_down(:,i,j),STD_CCA_down(:,i,j),'color',[0.5 0.5 0.5],'LineWidth',2);
    end
    if any(lstHrfAdd(:,1) == i)
        First_line = ['HRF, Corr: ' num2str(CORR_CCA(foo,j),'%0.2g')];
        Second_line =  [' MSE: ' num2str(MSE_CCA(foo,j),'%0.2g') ];
        xlabel({First_line;Second_line})
        foo = foo + 1;
        plot([0:1/fq:max(hrf.t_hrf+1/fq)],hrf.hrf_conc(:,j));
    end
    txt = ['ch ' num2str(i)];
    ylabel(txt);
    grid; ylim([ylim1 ylim2]);xlim([xlim1 xlim2]);
end
suptitle(['GLM with CCA - Subject # ' num2str(ss) ',  t_l_a_g= ' num2str(timelag) ' sec']) ;

