% plot a ref pair broke up by channel 

clear all;
addpath(genpath('D:\Data\Pipeline\npy\npy-matlab'))
day1 = 1;
day2 = 3;
clu_label1 = 50;
clu_label2 = 60;

fig_path = 'C:\Users\labadmin\Desktop\Neuron Tracking Pipeline\User version\Github Figure reproduce\F1'; %NEED CHANGE to local path
wave1 = readNPY(fullfile(fig_path,'ksproc_mean_waveforms1.npy'));
wave2 = readNPY(fullfile(fig_path,'ksproc_mean_waveforms2.npy'));
chann_map = readNPY(fullfile(fig_path,'channel_map.npy'))'+1;
chann_pos = load(fullfile(fig_path,'chan_pos.mat')).chann_pos;
clu = load(fullfile(fig_path,'clu.mat')).clu;
ishank = 1; 
clu_idx1 = find(clu(day1,ishank).clu == clu_label1); %cluster index 
clu_idx2 = find(clu(day2,ishank).clu == clu_label2);


peakWf1 = squeeze(wave1(clu_idx1,:,:)); %wf of one cluster
peakWf2 = squeeze(wave2(clu_idx2,:,:));


rowHalfRange = 2;
xMax = 40;
ylow = 0;
yhigh = 5000;
figName = 'test_spike';

plotWaves(peakWf1, peakWf2, chann_pos(:,1), chann_pos(:,2), 1,1, rowHalfRange, xMax, yhigh, ylow, figName);
legend('','Day 1 unit 50','Day 13 unit 60')
title(sprintf('AL032 shank 3 Day %d unit %d and Day 13 unit %d',day1,clu_label1,clu_label2))





%-----------------Helper function-----------------
function [yhigh, ylow] = plotWaves(peakWf1, peakWf2, xC, yC, rowFac, colFac, rowHalfRange, xMax, yhigh, ylow, figName)

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
hSep = colPos(2) - colPos(1); %x position step size
vSep = rowPos(2) - rowPos(1); %z position step size

% plot setting
colSpace = maxSample + tPadFrac*maxSample; % horizontal dimensions
xlow = -0.5*colSpace;
xhigh = (xMax/hSep)*colSpace*colFac;
% want 0:nt-1 to equal hSep
tscale = ((xhigh/xMax)*hSep)/colSpace; % points in hSep

tpts = (0:maxSample-1)*tscale; %time points

p_halfW = 0.5*max(tpts); %shift left right
p_halfH = 0.725*rowSpace; %shift up down
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
