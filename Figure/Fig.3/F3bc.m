% function emd_drift_figure

clear all

fig_path = 'C:\Users\labadmin\Desktop\Neuron Tracking Pipeline\User version\Github Figure reproduce\F3&4'; %NEED CHANGE to local path
load(fullfile(fig_path,'F3bc_data.mat'))

hasMatch = ~isnan(f2_emd_ind);
f1_matched = f1(hasMatch,:);
f2_matched = f2(f2_emd_ind(hasMatch),:);
diffZ = f2_matched(:,2) - f1_matched(:,2);
binsize = 4;
edges = (-100:binsize:100);
mode_diffZ = kernelModeEstimate(diffZ);
lim_diffZ = 15;  % 2*4, est err



%% line color by threshold
h2 = figure();
h2.Name = 'emd pairs';
h2.Units = 'centimeters';
h2.Position = [8,4,6,24];
% get(h2,'Position')

scatter(f1(:,1), f1(:,2), 50,'k','o'); hold on;
scatter(f2(:,1), f2(:,2),50,'k','o','filled'); hold on
for i = 1:length(f1_matched)
    u = f2_matched(i,1)-f1_matched(i,1);
    v = f2_matched(i,2)-f1_matched(i,2);

    % movement < mode + threshold
    if abs(diffZ(i)-mode_diffZ) < lim_diffZ
        if diffZ(i) < 0 %move down
            quiver(f1_matched(i,1),f1_matched(i,2),u,v,'Color','black','LineWidth', 1.5,'MaxHeadSize',0.5);
        else %move up = correct
            quiver(f1_matched(i,1),f1_matched(i,2),u,v,'Color','red','LineWidth', 1.5,'MaxHeadSize',0.5);
        end
    else % above threshold
        if diffZ(i) < 0
            quiver(f1_matched(i,1),f1_matched(i,2),u,v,'Color', 'black','LineStyle',':','LineWidth', 1,'MaxHeadSize',0.5);
        else
            quiver(f1_matched(i,1),f1_matched(i,2),u,v,'Color', 'red','LineStyle',':','LineWidth', 1,'MaxHeadSize',0.5);
        end
    end
    hold on;
end
xlim([-20,50]);
ylim([2890,3600]);
xticks([-20 0 20 40])
xtickangle(0)
ax = gca;
% set(ax,'YTick',[])
set(ax,'PlotBoxAspectRatio',[0.157,1,0.157])
set(ax,'TickDir','out');
set(ax,'box','on')
ax.FontSize = 16;
legend('Day 1 units','Day 2 units','shift up (short distance)','','','','','shift down (long distance)','','','','','','','','','','','','','','','shift up (long distance)')
xlabel('X (\mum)','FontSize',18,'FontWeight','Bold','Color','k')



% enlarged
figure()
scatter(f1(:,1), f1(:,2), 50,'k','o'); hold on;
scatter(f2(:,1), f2(:,2),50,'k','o','filled'); hold on
for i = 1:length(f1_matched)
    u = f2_matched(i,1)-f1_matched(i,1);
    v = f2_matched(i,2)-f1_matched(i,2);

    % movement < mode + threshold
    if abs(diffZ(i)-mode_diffZ) < lim_diffZ
        if diffZ(i) < 0 %move down
            quiver(f1_matched(i,1),f1_matched(i,2),u,v,'Color','black','LineWidth', 2,'MaxHeadSize',0.8);
        else %move up = correct
            quiver(f1_matched(i,1),f1_matched(i,2),u,v,'Color','red','LineWidth', 2,'MaxHeadSize',0.8);
        end
    else % above threshold
        if diffZ(i) < 0
            quiver(f1_matched(i,1),f1_matched(i,2),u,v,'Color', 'black','LineStyle',':','LineWidth', 1,'MaxHeadSize',0.8);
        else
            quiver(f1_matched(i,1),f1_matched(i,2),u,v,'Color', 'red','LineStyle',':','LineWidth', 1,'MaxHeadSize',0.8);
        end
    end
    hold on;
end
xlim([-20,40]);
ylim([3120,3225]);
ax = gca;
ax.FontSize = 16;
set(ax,'TickDir','out');
set(ax,'box','on')





%% Histogram
h3 = figure();
h3.Name = 'zdiff_his';
h3.Units = 'centimeters';
% h3.Position = [4,4,fwidth,fheight];
h3.Position = [20,10,15,10];

blue = [0.3 0.75 0.93];
purple = [0.72 0.27 1];
z_high = max(diffZ);
z_low = min(diffZ);
nbins = 60;
% edge = [z_low:4:z_high];
[hf,edge] = histcounts(diffZ,nbins);
r_upper = mode_diffZ + lim_diffZ;
r_lower = mode_diffZ - lim_diffZ;
hh = histfit(diffZ,nbins,'kernel'); hold on
% hh(1).FaceColor = blue; %histogram color
hh(1).FaceColor = 'flat';
hh(1).CData(1:21,:) = repmat([0 0 0],21,1); %out of threshold bars
hh(1).CData(22:27,:) = repmat([1 0 0],6,1); %in threshold bars
hh(1).CData(28:nbins,:) = repmat([0 0 0],33,1);
hh(1).EdgeColor = 'none'; %remove histogram edge
hh(2).Color = blue; %kernel fit color
xline(mode_diffZ,'--','color',blue,'LineWidth',2.5); hold on
xline(12,'--','color','b','LineWidth',2.5); hold on
xline(mode_diffZ+lim_diffZ,'r:','LineWidth',1.5); hold on
xline(mode_diffZ-lim_diffZ,'r:','LineWidth',1.5);

ax = gca;
ax.FontSize = 18;
set(ax,'box','off')
set(ax,'TickDir','out'); %tickmark towards outside
xlim([-50,60])
xticks([-50:25:50,60])
legend('implausible distance','acceptable distance')
xlabel('EMD matched Vertical-shift (\mum)','FontSize',18,'FontWeight','bold','Color','k')
ylabel('Number of matches','FontSize',18,'FontWeight','bold','Color','k')



