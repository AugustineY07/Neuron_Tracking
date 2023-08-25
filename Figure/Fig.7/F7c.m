% Plot for poster figure 2, example RF

clear all

fig_path = 'C:\Users\labadmin\Desktop\Neuron Tracking Pipeline\User version\Github Figure reproduce\F2'; %NEED CHANGE to local path
load(fullfile(fig_path,'F2c_data.mat'))


idx_unit1 = 1; %first pair
idx_unit2 = 22; %second pair
clu1 = ref12(idx_unit1,6);
clu2 = ref12(idx_unit1,5);
clu3 = ref15(idx_unit2,6);
clu4 = ref15(idx_unit2,5);
rsp1 = rsp_sig{1,1}(ref12(idx_unit1,2),:); %rsp for cluster
rsp2 = rsp_sig{2,1}(ref12(idx_unit1,1),:);
rsp3 = rsp_sig{1,1}(ref15(idx_unit2,2),:); 
rsp4 = rsp_sig{5,1}(ref15(idx_unit2,1),:); 

psth1 = PSTH_sig{1,1}(ref12(idx_unit1,2),:); %psth for cluster
psth2 = PSTH_sig{2,1}(ref12(idx_unit1,1),:);
psth3 = PSTH_sig{1,1}(ref15(idx_unit2,2),:); 
psth4 = PSTH_sig{5,1}(ref15(idx_unit2,1),:); 


Dblue = [0 0.45 0.74];
Dred = [0.64 0.08 0.18];

figure()
subplot(2,2,1)
[hAx,hLine1,hLine2] = plotyy(1:numel(rsp1),rsp1, 1:numel(rsp2),rsp2);
set(hAx,'FontSize', 30); %tick font
set(hAx,'Box', 'off'); %remove tick box
set(hAx,'TickDir','out'); %tickmark towards outside
p1 = get(hAx,'position');
set(hAx,'position',[0.155 0.5838 0.3347 0.3412])
hLine1.LineWidth = 2;
hLine1.Color = Dblue;
hLine2.LineWidth = 2;
hLine2.LineStyle = '-.';
hLine2.Color = Dred;
set(hAx(1),'YTick',[0:10:50])
set(hAx(1),'YColor',Dblue);
set(hAx(2),'YTick',[])
set(hAx(2),'YColor',Dred);
xlabel('Stimulus id','FontSize',30,'FontWeight','Bold','Color','k')
ylabel('Spike rate (day 1)','FontSize',30,'FontWeight','Bold','Color',Dblue) % left y-axis
legend({sprintf('day 1, unit %d',clu1),sprintf('day 2, unit %d',clu2)},'Location','northeast','Orientation','vertical','FontSize',30)
h1 = subtitle(sprintf('corr\\_vfp = %.2f', sim_rsp_sig{2,1}(ref12(1,2), ref12(1,1)))); 
h1.FontSize = 26;
pbaspect([1 0.7 0.7])

subplot(2,2,3)
[hAx2,hLine3,hLine4] = plotyy(1:numel(psth1),psth1, 1:numel(psth2),psth2);
set(hAx2,'FontSize', 30); %tick font
set(hAx2,'Box', 'off'); %remove tick box
set(hAx2,'TickDir','out'); %tickmark towards outside
p2 = get(hAx2,'position');
set(hAx2,'position',[0.155 0.11 0.3347 0.3412])
hLine3.LineWidth = 2;
hLine3.Color = Dblue;
hLine4.LineWidth = 2;
hLine4.LineStyle = '-.';
hLine4.Color = Dred ;
set(hAx2,'XTick',[0:500:2000])
set(hAx2(1),'YTick',[0:0.5:1.5])
set(hAx2(1),'YColor',Dblue);
set(hAx2(2),'YTick',[])
set(hAx2(2),'YColor',Dred);
xlabel('Time (ms)','FontSize',30,'FontWeight','Bold','Color','k')
ylabel('Normalized spike rate (day 1)','FontSize',30,'FontWeight','Bold','Color',Dblue) % left y-axis
legend({sprintf('day 1, unit %d',clu1),sprintf('day 2, unit %d',clu2)},'Location','northeast','Orientation','vertical','FontSize',30)
h2 = subtitle(sprintf('corr\\_PSTH = %.2f', sim_PSTH_sig{2,1}(ref12(1,2), ref12(1,1))));
h2.FontSize = 26;

subplot(2,2,2)
[hAx3,hLine5,hLine6] = plotyy(1:numel(rsp3),rsp3, 1:numel(rsp4),rsp4);
set(hAx3,'FontSize', 30); %tick font
set(hAx3,'Box', 'off'); %remove tick box
set(hAx3,'TickDir','out'); %tickmark towards outside
p3 = get(hAx3,'position');
set(hAx3,'position',[0.545 0.5838 0.3347 0.3412])
hLine5.LineWidth = 2;
hLine5.Color = Dblue;
hLine6.LineWidth = 2;
hLine6.LineStyle = '-.';
hLine6.Color = Dred;
set(hAx3(1),'YLim',[0 50])
set(hAx3(1),'YTick',[0:10:50])
set(hAx3(1),'YColor',Dblue);
set(hAx3(2),'YLim',[0 50])
set(hAx3(2),'YTick',[0:10:50])
set(hAx3(2),'YColor',Dred);
xlabel('Stimulus id','FontSize',30,'FontWeight','Bold','Color','k')
ylabel(hAx3(1),'Spike rate (day 1)','FontSize',30,'FontWeight','Bold','Color',Dblue) % left y-axis
ylabel(hAx3(2),'Spike rate (day 48)','FontSize',30,'FontWeight','Bold','Color',Dred) % right y-axis
legend({sprintf('day 1, unit %d',clu3),sprintf('day 48, unit %d',clu4)},'Location','northeast','Orientation','vertical','FontSize',30)
h3 = subtitle(sprintf('corr\\_vfp = %.2f', sim_rsp_sig{5,1}(ref15(idx_unit2,2), ref15(idx_unit2,1)))); %'fontsize',16
h3.FontSize = 26;

subplot(2,2,4)
[hAx4,hLine7,hLine8] = plotyy(1:numel(psth3),psth3, 1:numel(psth4),psth4);
set(hAx4,'FontSize', 30); %tick font
set(hAx4,'Box', 'off'); %remove tick box
set(hAx4,'TickDir','out'); %tickmark towards outside
p4 = get(hAx4,'position');
set(hAx4,'position',[0.545 0.11 0.3347 0.3412])
hLine7.LineWidth = 2;
hLine7.Color = Dblue;
hLine8.LineWidth = 2;
hLine8.LineStyle = '-.';
hLine8.Color = Dred;
set(hAx4,'XTick',[0:500:2000])
set(hAx4(1),'YLim',[0 1.5])
set(hAx4(1),'YTick',[0:0.5:1.5])
set(hAx4(1),'YColor',Dblue);
set(hAx4(2),'YLim',[0 1.5])
set(hAx4(2),'YTick',[0:0.5:1.5])
set(hAx4(2),'YColor',Dred);
xlabel('Time (ms)','FontSize',30,'FontWeight','Bold','Color','k')
ylabel(hAx4(1),'Normalized spike rate (day 1)','FontSize',30,'FontWeight','Bold','Color',Dblue) % left y-axis
ylabel(hAx4(2),'Normalized spike rate (day 48)','FontSize',30,'FontWeight','Bold','Color',Dred) % right y-axis
legend({sprintf('day 1, unit %d',clu3),sprintf('day 48, unit %d',clu4)},'Location','northeast','Orientation','vertical','FontSize',30)
h4 = subtitle(sprintf('corr\\_PSTH = %.2f', sim_PSTH_sig{5,1}(ref15(idx_unit2,2), ref15(idx_unit2,1))));
h4.FontSize = 26;

