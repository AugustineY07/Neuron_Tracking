function plot_waveform(input, clu_label1, clu_label2)
% plot waveform traces of matched units

chan_pos = input.chan_pos;
chan_map = input.chan_map;
input_path = input.input_path;
data_path1 = input.data_path1;
data_path2 = input.data_path2;
z_step = input.zStep;
x_step = input.xStep;
wave1 = readNPY(fullfile(input_path,data_path1,input.wf_name));
wave2 = readNPY(fullfile(input_path,data_path2,input.wf_name));

% if the number of channels in mw differ from chan_map, assume mw is the
% original data including all channels, and select only those included in
% the sort.
[nChanPos, ~] = size(chan_pos);
[~, nChanMW, ~] = size(wave1);

if nChanPos < nChanMW
    wave1 = wave1(:,chan_map+1,:);
    wave2 = wave2(:,chan_map+1,:);
end


%read channel info
xC = chan_pos(:,1); yC = chan_pos(:,2);

peakWf1 = squeeze(wave1(clu_label1,:,:)); %wf of one cluster, get by the cluster label
peakWf2 = squeeze(wave2(clu_label2,:,:));

rowHalfRange = 4; %how many rows to plot
xMax = 40;
ylow = 0;
yhigh = 5000;
figName = 'example_spike';

plotWaves(peakWf1, peakWf2, xC, yC, x_step, z_step, 1,1, rowHalfRange, xMax, yhigh, ylow, figName);
legend('',sprintf('Unit %d',clu_label1),sprintf('Unit %d',clu_label2))
% title(sprintf('AL032 shank 3 Day %d unit %d and Day 13 unit %d',day1,clu_label1,clu_label2))
end


% rowFac=1; colFac=1;
function [yhigh, ylow] = plotWaves(peakWf1, peakWf2, xC, yC,hSep,vSep, rowFac, colFac, rowHalfRange, xMax, yhigh, ylow, figName)

% plot waveforms in rough coordinate positions
bAxis = 0;  % whether or not to include axes in the figure

% figure dimensions
fwidth = 10; % figure width, in cm
fheight = 1.5*fwidth; % adjust to set aspect ratio
rowSpace = 150; % uV
tPadFrac = 0.3; % increase to increase the space between columns



[nChan,nt] = size(peakWf1);
maxSample = nt;
colPos = unique(xC); %all x positions
rowPos = unique(yC); %all y positions

% plot setting
colSpace = maxSample + tPadFrac*maxSample; % horizontal dimensions
xlow = -0.5*colSpace;
xhigh = (xMax/hSep)*colSpace*colFac;
% want 0:nt-1 to equal hSep
tscale = ((xhigh/xMax)*hSep)/colSpace; % points in hSep

tpts = (0:maxSample-1)*tscale; %time points

p_halfW = 0.5*max(tpts); %shift background box left right
p_halfH = 0.9*rowSpace; %shift background box up down
h = figure('Name',figName,'Units','centimeters','Position',[4,4,fwidth,fheight]);


for i = 1:nChan
    currCol = find(colPos==xC(i)) - 1;
    currRow = find(rowPos==yC(i)) - 1;
    currWave1 = squeeze(peakWf1(i,:));
    currWave2 = squeeze(peakWf2(i,:));
    
    currX = tpts + currCol*colFac*colSpace;
    currY = currWave1 + currRow*rowFac*rowSpace;
    currY2 = currWave2 + currRow*rowFac*rowSpace;

    currX = currX(1:maxSample);
    currY = currY(1:maxSample);
    currY2 = currY2(1:maxSample);

    xOff = currCol*colFac*colSpace + p_halfW;
    yOff = currRow*rowFac*rowSpace;
    pgon = polyshape([xOff-p_halfW  xOff+p_halfW  xOff+p_halfW  xOff-p_halfW], ...
                                [yOff-p_halfH  yOff-p_halfH yOff+p_halfH yOff+p_halfH]);
    plot(pgon,'FaceColor','white','FaceAlpha',1,'EdgeColor','Red','EdgeAlpha',0.1); hold on
    plot(currX, currY,'Color',[0 0.45 0.74],'LineWidth',1.5); hold on
    plot(currX, currY2,'-.','Color',[0.64 0.08 0.18],'LineWidth',1.5);
    
    hold on;
end


% first term is offset to middle of trace in column 0
xtick_arr = 0.5*maxSample*tscale + (0:numel(colPos)-1)*colFac*colSpace;
xticks(xtick_arr);

xticklabels((0:numel(xticks)-1)*hSep);
% y ticks centered on traces
ytick_arr = ((0:numel(rowPos)-1))*rowFac*rowSpace;
if vSep > 6
    yticks(ytick_arr);
    yticklabels((0:numel(yticks)-1)*vSep);
else
    ytick_arr = ytick_arr(1:3:end);
    yticks(ytick_arr);
    yticklabels((0:numel(yticks)-1)*3*vSep);
end

xlim([-0.25*colSpace, (numel(colPos)-1)*colSpace*colFac + 1.25*colSpace])

if rowHalfRange > -1
    % calculate the y limits
    minV_allSite = squeeze(min(peakWf1,[],2));
    [~, maxInd] = min(minV_allSite);
    maxRow = floor((yC(maxInd)-min(yC))/vSep); % usually, yC will be even multiples of vSep
    ylow = -0.5*rowSpace + (maxRow-rowHalfRange)*rowSpace*rowFac;
    yhigh = (maxRow+rowHalfRange)*rowSpace*rowFac + 0.5*rowSpace;
    ylim([ylow,yhigh]);
elseif (yhigh ~= 0)
    % use the input values
    ylim([ylow,yhigh]);
else
    % use whole range of sites in waves
    ylim([-0.5*rowSpace, (numel(rowPos)-1)*rowSpace*rowFac + 1.5*rowSpace])
end


if ~bAxis
    axis off
end
hold off;

end

%%
% % figure dimensions
% fwidth = 10; % figure width, in cm
% fheight = 1.5*fwidth; % adjust to set aspect ratio
% rowSpace = 600; % uV, box height
% tPadFrac = 0.3; % increase to increase the space between columns
% 
% % plot setting
% [nChan,maxSample] = size(peakWf1);
% colPos = unique(xC); %all x positions
% rowPos = unique(yC); %all y positions
% hSep = colPos(2) - colPos(1); %x position step size
% vSep = rowPos(2) - rowPos(1); %z position step size
% 
% %cluster 1
% offset = max(peakWf1,[],2) - min(peakWf1,[],2);
% [~, maxInd] = max(offset);
% currentRow = (yC(maxInd)-min(yC))/vSep; %current row, count from 0
% startRow = max([currentRow-2,0]);
% startSite = find(yC == startRow*vSep+min(yC),1,'first');
% endRow = min([currentRow+2,(max(yC)-min(yC))/vSep]); %number of rows below
% endSite = find(yC == endRow*vSep+min(yC),1,'last');
% peakRow = startRow; %relative location of peak row, use for alignment
% numSite = (endRow - startRow + 1)*2; %number of sites plotted, might involve 1 reference site
% x_site = repmat([min(xC),max(xC)],1,5);
% rows = repelem([startRow:1:endRow],2);
% siteCount = 0;
% for isite = 1:numSite %find relevant site waveforms
%     this_site = find(xC == x_site(isite) & yC == rows(isite)*vSep+min(yC));
%     if ~isempty(this_site)
%         siteCount = siteCount + 1;
%         site(siteCount) = this_site;
%         plotSite(siteCount,:) = peakWf1(site(siteCount),:); %peak sites of this cluster
%         xC_sub1 = xC(site); %x locations
%         yC_sub1 = yC(site); %y locations
%     end
% end
% 
% %cluster 2
% offset = max(peakWf2,[],2) - min(peakWf2,[],2);
% [~, maxInd] = max(offset);
% currentRow = (yC(maxInd)-min(yC))/vSep; %current row, count from 0
% startRow = max([currentRow-2,0]);
% startSite = find(yC == startRow*vSep+min(yC),1,'first');
% endRow = min([currentRow+2,(max(yC)-min(yC))/vSep]); %number of rows below
% endSite = find(yC == endRow*vSep+min(yC),1,'last');
% peakRow = startRow; %relative location of peak row, use for alignment
% 
% numSite = (endRow - startRow + 1)*2;
% x_site = repmat([min(xC),max(xC)],1,5);
% rows = repelem([startRow:1:endRow],2);
% siteCount = 0;
% for isite = 1:numSite %find relevant site index
%     this_site = find(xC == x_site(isite) & yC == rows(isite)*vSep+min(yC));
%     if ~isempty(this_site)
%         siteCount = siteCount + 1;
%         site(siteCount) = this_site;
%         plotSite(siteCount,:) = peakWf2(site(siteCount),:); %peak sites of this cluster
%         xC_sub2 = xC(site); %x locations
%         yC_sub2 = yC(site); %y locations
%     end
% end
% 
% 
% % plot setting
% xMax = 40;
% rowFac = 1; colFac = 1;
% colSpace = maxSample + tPadFrac*maxSample; % start of 2nd box in row
% xlow = -0.5*colSpace;
% xhigh = (xMax/hSep)*colSpace*colFac;
% % want 0:nt-1 to equal hSep
% tscale = ((xhigh/xMax)*hSep)/colSpace; % points in hSep = 1
% tpts = (0:maxSample-1)*tscale; %time points
% p_halfW = 0.5*max(tpts); %shift left right, x center of box 1
% p_halfH = 0.5*rowSpace; %shift up down, y center of box 1
% 
% 
% 
% figure()
% common_row = min(size(plotSite,1),size(plotSite,1));
% for i = 1:common_row %plot common rows
%     currCol = find(colPos==xC_sub1(i)) - 1; %current box row and col
%     currRow = find(rowPos==yC_sub1(i)) - 1;
%     currWave1 = plotSite(i,:); %wf on a single channel
%     currWave2 = plotSite(i,:);
% 
%     currX = tpts + currCol*colFac*colSpace; %all time point x location
%     currY = currWave1 + currRow*rowFac*rowSpace; %all wf amp y location
%     currY2 = currWave2 + currRow*rowFac*rowSpace;
% 
%     currX = currX(1:maxSample);
%     currY = currY(1:maxSample);
%     currY2 = currY2(1:maxSample);
% 
%     xOff = currCol*colFac*colSpace + p_halfW;
%     yOff = currRow*rowFac*rowSpace;
%     pgon = polyshape([xOff-p_halfW  xOff+p_halfW  xOff+p_halfW  xOff-p_halfW], ...
%         [yOff-p_halfH  yOff-p_halfH yOff+p_halfH yOff+p_halfH]);
%     plot(pgon,'FaceColor','white','FaceAlpha',1,'EdgeColor','Red','EdgeAlpha',0.1); hold on
%     plot(currX, currY,'Color',[0 0.45 0.74],'LineWidth',1.5); hold on
%     plot(currX, currY2,'-.','Color',[0.64 0.08 0.18],'LineWidth',1.5); hold on;
% end
% % plot extra rows
% for i = common_row+1:max(size(plotSite,1),size(plotSite,1))
%     if i <= size(plotSite,1)
%         currCol = find(colPos==xC_sub1(i)) - 1; %current box row and col
%         currRow = find(rowPos==yC_sub1(i)) - 1;
%         currWave1 = plotSite(i,:); %wf on a single channel
% 
%         currX = tpts + currCol*colFac*colSpace; %all time point x location
%         currY = currWave1 + currRow*rowFac*rowSpace; %all wf amp y location
% 
%         currX = currX(1:maxSample);
%         currY = currY(1:maxSample);
% 
%         xOff = currCol*colFac*colSpace + p_halfW;
%         yOff = currRow*rowFac*rowSpace;
%         pgon = polyshape([xOff-p_halfW  xOff+p_halfW  xOff+p_halfW  xOff-p_halfW], ...
%             [yOff-p_halfH  yOff-p_halfH yOff+p_halfH yOff+p_halfH]);
%         plot(pgon,'FaceColor','white','FaceAlpha',1,'EdgeColor','Red','EdgeAlpha',0.1); hold on
%         plot(currX, currY,'Color',[0 0.45 0.74],'LineWidth',1.5); hold on
%     elseif i <= size(plotSite,1)
%         currCol = find(colPos==xC_sub2(i)) - 1; %current box row and col
%         currRow = find(rowPos==yC_sub2(i)) - 1;
%         currWave2 = plotSite(i,:);
% 
%         currX = tpts + currCol*colFac*colSpace; %all time point x location
%         currY2 = currWave2 + currRow*rowFac*rowSpace;
% 
%         currX = currX(1:maxSample);
%         currY2 = currY2(1:maxSample);
% 
%         xOff = currCol*colFac*colSpace + p_halfW;
%         yOff = currRow*rowFac*rowSpace;
%         pgon = polyshape([xOff-p_halfW  xOff+p_halfW  xOff+p_halfW  xOff-p_halfW], ...
%             [yOff-p_halfH  yOff-p_halfH yOff+p_halfH yOff+p_halfH]);
%         plot(pgon,'FaceColor','white','FaceAlpha',1,'EdgeColor','Red','EdgeAlpha',0.1); hold on
%         plot(currX, currY2,'-.','Color',[0.64 0.08 0.18],'LineWidth',1.5); hold on;
%     end
% end
% 
% legend('',sprintf('Unit %d',clu_label1),sprintf('Unit %d',clu_label2))
% ax = gca;
% ax.FontSize = 16; %tick font
% ax.Box = 'off'; %remove tick box
% set(ax,'TickDir','out'); %tickmark towards outside
% 
% axis off

