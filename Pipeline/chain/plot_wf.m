function L2_weight = plot_wf(all_input, full_chain, L2_weight, chan_pos, numData, ichain, id)
% figure dimensions
fwidth = 10; % figure width, in cm
fheight = 1.5*fwidth; % adjust to set aspect ratio
rowSpace = 600; % uV, box height
tPadFrac = 0.3; % increase to increase the space between columns


input_path = all_input(id).input.input_path;
data_path1 = all_input(id).input.data_path1;
data_path2 = all_input(id).input.data_path2;
wave1 = readNPY(fullfile(input_path,data_path1,'ksproc_mean_waveforms.npy'));
wave2 = readNPY(fullfile(input_path,data_path2,'ksproc_mean_waveforms.npy'));
%read channel info
xC = chan_pos(:,1); yC = chan_pos(:,2);

%current clusters on the chain
clu_label1 = full_chain(ichain,id);
clu_label2 = full_chain(ichain,id+1);
peakWf1 = squeeze(wave1(clu_label1,:,:)); %wf of one cluster, get by the cluster label
peakWf2 = squeeze(wave2(clu_label2,:,:));

% plot setting
[nChan,maxSample] = size(peakWf1);
colPos = unique(xC); %all x positions
rowPos = unique(yC); %all y positions
hSep = colPos(2) - colPos(1); %x position step size = 32
vSep = rowPos(2) - rowPos(1); %z position step size = 15

%cluster 1
offset = max(peakWf1,[],2) - min(peakWf1,[],2);
[~, maxInd] = max(offset);
currentRow = (yC(maxInd)-min(yC))/vSep; %current row, count from 0
startRow = max([currentRow-2,0]);
startSite = find(yC == startRow*vSep+min(yC),1,'first');
endRow = min([currentRow+2,(max(yC)-min(yC))/vSep]); %number of rows below
endSite = find(yC == endRow*vSep+min(yC),1,'last');
peakRow(id) = startRow; %relative location of peak row, use for alignment
numSite = (endRow - startRow + 1)*2; %number of sites plotted, might involve 1 reference site
x_site = repmat([min(xC),max(xC)],1,5);
rows = repelem([startRow:1:endRow],2);
siteCount = 0;
for isite = 1:numSite %find relevant site waveforms
    this_site = find(xC == x_site(isite) & yC == rows(isite)*vSep+min(yC));
    if ~isempty(this_site)
        siteCount = siteCount + 1;
        site(siteCount) = this_site;
        plotSite{id}(siteCount,:) = peakWf1(site(siteCount),:); %peak sites of this cluster
        xC_sub1 = xC(site); %x locations
        yC_sub1 = yC(site); %y locations
    end
end

%cluster 2
offset = max(peakWf2,[],2) - min(peakWf2,[],2);
[~, maxInd] = max(offset);
currentRow = (yC(maxInd)-min(yC))/vSep; %current row, count from 0
startRow = max([currentRow-2,0]);
startSite = find(yC == startRow*vSep+min(yC),1,'first');
endRow = min([currentRow+2,(max(yC)-min(yC))/vSep]); %number of rows below
endSite = find(yC == endRow*vSep+min(yC),1,'last');
peakRow(id+1) = startRow; %relative location of peak row, use for alignment

numSite = (endRow - startRow + 1)*2;
x_site = repmat([min(xC),max(xC)],1,5);
rows = repelem([startRow:1:endRow],2);
siteCount = 0;
for isite = 1:numSite %find relevant site index
    this_site = find(xC == x_site(isite) & yC == rows(isite)*vSep+min(yC));
    if ~isempty(this_site)
        siteCount = siteCount + 1;
        site(siteCount) = this_site;
        plotSite{id+1}(siteCount,:) = peakWf2(site(siteCount),:); %peak sites of this cluster
        xC_sub2 = xC(site); %x locations
        yC_sub2 = yC(site); %y locations
    end
end


% plot setting
xMax = 40;
rowFac = 1; colFac = 1;
colSpace = maxSample + tPadFrac*maxSample; % start of 2nd box in row
xlow = -0.5*colSpace;
xhigh = (xMax/hSep)*colSpace*colFac;
% want 0:nt-1 to equal hSep
tscale = ((xhigh/xMax)*hSep)/colSpace; % points in hSep = 1
tpts = (0:maxSample-1)*tscale; %time points
p_halfW = 0.5*max(tpts); %shift left right, x center of box 1
p_halfH = 0.5*rowSpace; %shift up down, y center of box 1


%plot wfs
subplot(1,numData-1,id)
common_row = min(size(plotSite{id},1),size(plotSite{id+1},1));
for i = 1:common_row %plot common rows
    currCol = find(colPos==xC_sub1(i)) - 1; %current box row and col
    currRow = find(rowPos==yC_sub1(i)) - 1;
    currWave1 = plotSite{id}(i,:); %wf on a single channel
    currWave2 = plotSite{id+1}(i,:);

    currX = tpts + currCol*colFac*colSpace; %all time point x location
    currY = currWave1 + currRow*rowFac*rowSpace; %all wf amp y location
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
    plot(currX, currY2,'-.','Color',[0.64 0.08 0.18],'LineWidth',1.5); hold on;
end
% plot extra rows
for i = common_row+1:max(size(plotSite{id},1),size(plotSite{id+1},1))
    if i <= size(plotSite{id},1)
        currCol = find(colPos==xC_sub1(i)) - 1; %current box row and col
        currRow = find(rowPos==yC_sub1(i)) - 1;
        currWave1 = plotSite{id}(i,:); %wf on a single channel

        currX = tpts + currCol*colFac*colSpace; %all time point x location
        currY = currWave1 + currRow*rowFac*rowSpace; %all wf amp y location

        currX = currX(1:maxSample);
        currY = currY(1:maxSample);

        xOff = currCol*colFac*colSpace + p_halfW;
        yOff = currRow*rowFac*rowSpace;
        pgon = polyshape([xOff-p_halfW  xOff+p_halfW  xOff+p_halfW  xOff-p_halfW], ...
            [yOff-p_halfH  yOff-p_halfH yOff+p_halfH yOff+p_halfH]);
        plot(pgon,'FaceColor','white','FaceAlpha',1,'EdgeColor','Red','EdgeAlpha',0.1); hold on
        plot(currX, currY,'Color',[0 0.45 0.74],'LineWidth',1.5); hold on
    elseif i <= size(plotSite{id+1},1)
        currCol = find(colPos==xC_sub2(i)) - 1; %current box row and col
        currRow = find(rowPos==yC_sub2(i)) - 1;
        currWave2 = plotSite{id+1}(i,:);

        currX = tpts + currCol*colFac*colSpace; %all time point x location
        currY2 = currWave2 + currRow*rowFac*rowSpace;

        currX = currX(1:maxSample);
        currY2 = currY2(1:maxSample);

        xOff = currCol*colFac*colSpace + p_halfW;
        yOff = currRow*rowFac*rowSpace;
        pgon = polyshape([xOff-p_halfW  xOff+p_halfW  xOff+p_halfW  xOff-p_halfW], ...
            [yOff-p_halfH  yOff-p_halfH yOff+p_halfH yOff+p_halfH]);
        plot(pgon,'FaceColor','white','FaceAlpha',1,'EdgeColor','Red','EdgeAlpha',0.1); hold on
        plot(currX, currY2,'-.','Color',[0.64 0.08 0.18],'LineWidth',1.5); hold on;
    end
end


legend('',sprintf('Day %d unit %d',id,clu_label1),sprintf('Day %d unit %d',id+1,clu_label2))
subtitle({sprintf('Chain %d: Dataset %d u%d and %d u%d',ichain,id,clu_label1,id+1,clu_label2), sprintf('L2 = %.2f',L2_weight(ichain,id))},'FontSize',18,'FontWeight','Bold','Color','k')
ax = gca;
ax.FontSize = 16; %tick font
ax.Box = 'off'; %remove tick box
set(ax,'TickDir','out'); %tickmark towards outside

axis off


end