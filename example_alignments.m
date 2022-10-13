
close all
clear all

k = 1;
subject = 'AL032'; %one mouse recording folder 
date = '2020-03-31'; %recording form a particular day
myKsDirSorting = '\\znas.cortexlab.net\Subjects\AL032\2020-03-31\ephys_K1\sorting'; %sorting folder
probe_number = 1; 



%% Define folders

subjectsFolder = 'E:\'; %subject = a mouse 

alignDir = fullfile(subjectsFolder, subject, date, 'alignments'); %define aligned folder
%driftPlotFolder = 'C:\Users\annaL\Documents\PhD\ephys_results\driftPlots';
%probeDataPlotFolder = 'C:\Users\annaL\Documents\PhD\ephys_results\probeDataPlots';

mkdir(alignDir)


%% Get basic info

[tags, hasEphys] = getEphysTags(subject, date);
% determine what exp nums exist
[expNums, blocks, hasBlock, pars, isMpep, tl, hasTimeline] = ...
    dat.whichExpNums(subject, date);
TLexp = expNums(hasTimeline);
TLexp = TLexp(end);
useFlipper = true;

%% align times (timeline to ephys)

% for any ephys, load the sync data
if hasEphys
    for t = 1:length(tags)
        if isempty(tags{t})
            [~, pdFlips, allET] = loadSyncChronic(subject, date);
        else
            [~, pdFlips, allET] = loadSyncChronic(subject, date, tags{t});
        end
        if useFlipper
         %   ephysFlips{t} = allET{7}{1};
         ephysFlips{t} = allET;%{7}{1};
        else
            ephysFlips{t} = pdFlips;
        end
    end
end

% synchronize multiple ephys to each other
if hasEphys
    if length(tags)>1
        for t2 = 2:length(tags)
            fprintf(1, 'correct ephys %s to %s\n', tags{t2}, tags{1});
            [~, b] = makeCorrection(ephysFlips{1}, ephysFlips{t2}, false);
            writeNPY(b, fullfile(alignDir, sprintf('correct_ephys_%s_to_ephys_%s.npy', tags{t2}, tags{1})));
        end
    end
end

% detect sync events from timelines
tlFlips = {};
for e = 1:length(expNums)
    if hasTimeline(e)
        Timeline = tl{e};
        tt = Timeline.rawDAQTimestamps;
        if useFlipper
            evTrace = Timeline.rawDAQData(:, strcmp({Timeline.hw.inputs.name}, 'flipper'));
            evT = schmittTimes(tt, evTrace, [3 4]); % all flips, both up and down
        else
            evTrace = Timeline.rawDAQData(:, strcmp({Timeline.hw.inputs.name}, 'photoDiode'));
            evT = schmittTimes(tt, evTrace, [3 4]); % all flips, both up and down
            evT = evT([true; diff(evT)>0.2]);
        end
        tlFlips{e} = evT;
    end
end

% match up ephys and timeline events:
% algorithm here is to go through each timeline available, figure out
% whether the events in timeline align with any of those in the ephys. If
% so, we have a conversion of events in that timeline into ephys
%
% Only align to the first ephys recording, since the other ones are aligned
% to that 
if hasEphys
    ef = ephysFlips{1};
    if useFlipper && ef(1)<0.001
        % this happens when the flipper was in the high state to begin with
        % - a turning on event is registered at the first sample. But here
        % we need to drop it. 
        ef = ef(2:end);
    end
    for e = 1:length(expNums)
        if hasTimeline(e)
            fprintf('trying to correct timeline %d to ephys\n', expNums(e));

            tlT = tlFlips{e};

            success=false;
            if length(tlT)==length(ef)
                % easy case: the two are exactly coextensive
                [~,b] = makeCorrection(ef, tlT, true);
                success = true;
            elseif length(tlT)<length(ef) && ~isempty(tlT)
                [~,b,success] = findCorrection(ef, tlT, false);
            elseif length(tlT)>length(ef) && ~isempty(tlT)
                [~,a,success] = findCorrection(tlT, ef, false);
                b = [1/a(1); -a(2)/a(1)];
            end
            if success
                writeNPY(b, fullfile(alignDir, ...
                    sprintf('correct_timeline_%d_to_ephys_%s.npy', ...
                    expNums(e), tags{1})));
                fprintf('success\n');
            else
                fprintf('could not correct timeline to ephys\n');
            end
        end
    end
end
%%
% match up blocks and mpeps to timeline in order:
% want to connect each block or mpep with part of a timeline. So go through
% each of these in order, looking through the timelines in sequence (of
% what hasn't already been matched) looking for a match. 
lastTimes = zeros(1,length(expNums));
for e = 1:length(expNums)
    if hasBlock(e)
        for eTL = 1:length(expNums)
            if hasTimeline(eTL)
                fprintf('trying to correct block %d to timeline %d\n', expNums(e), expNums(eTL));
                if useFlipper
                    % didn't get photodiode flips above, so get them now
                    Timeline = tl{eTL};
                    tt = Timeline.rawDAQTimestamps;
                    evTrace = Timeline.rawDAQData(:, strcmp({Timeline.hw.inputs.name}, 'photoDiode'));
                    pdT = schmittTimes(tt, evTrace, [1.5 2.8]);
                else
                    pdT = tlFlips{eTL};
                end
                block = blocks{e};
                sw = block.stimWindowUpdateTimes; 
                sw = sw(2:end); % sometimes need this? Why? how did sw
                % get an extra event at the beginning? 
                
                success = false;
                if length(sw)<=length(pdT) && length(sw)>1
                    [~,b,success,actualTimes] = findCorrection(pdT, sw, false);
                end
                if success
                    writeNPY(b, fullfile(alignDir, ...
                        sprintf('correct_block_%d_to_timeline_%d.npy', ...
                        expNums(e), expNums(eTL))));
                    writeNPY(actualTimes, fullfile(alignDir, ...
                        sprintf('block_%d_sw_in_timeline_%d.npy', ...
                        expNums(e), expNums(eTL))));
                    fprintf('  success\n');
                    lastTimes(eTL) = actualTimes(end);
                else
                    fprintf('  could not correct block %d to timeline %d\n', expNums(e), expNums(eTL));
                end
            end
        end
    
    end
end

%% get visual responses
f_neg = block.events.stimulusOnTimes;
first_exp=expNums(find(hasBlock));
first_exp = num2str(first_exp(2));
al = readNPY(fullfile(alignDir, sprintf('correct_block_%s_to_timeline_%d.npy', first_exp, TLexp)));
bTLtoMaster = readNPY(fullfile(alignDir, sprintf('correct_timeline_%d_to_ephys_%s.npy', TLexp, tags{probe_number})));
eventTimes = f_neg;
eventTimes = applyCorrection(eventTimes, al);
eventTimes = applyCorrection(eventTimes, bTLtoMaster);
[spikeTimes, spikeAmps, spikeDepths, spikeSites] = ksDriftmap(myKsDirSorting);
depthBinSize = 80; % in units of the channel coordinates, in this case µm
timeBinSize = 0.01; % seconds
bslWin = [-0.2 0]; % window in which to compute "baseline" rates for normalization
psthType = 'norm'; % show the normalized version
eventName = 'stimulus onset'; % for figure labeling
window = [-1 3];
[timeBins, depthBins, allP, ~] = psthByDepth(spikeTimes, spikeDepths, ...
    depthBinSize, timeBinSize, eventTimes, window, bslWin);

figure;
plotPSTHbyDepth(timeBins, depthBins, allP, eventName, psthType);