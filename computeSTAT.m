% find the preprocessed data binary file
fprintf('Loading preprocessed data and preparing output file... \n');
rootP = 'D:\allKsOutput\Science\A3\Preprocessed'; % the preprocessed data binary file is in this folder
fs          = dir(fullfile(rootP, '*.bin'));
fsfile = fullfile(rootP, fs.name);
% find channel positions in the shank 
Channel_x = h5read('D:\allKsOutput\Science\science.h5','/Dataset3/Channel_x');
Channel_y = h5read('D:\allKsOutput\Science\science.h5','/Dataset3/Channel_y');
channelPosition = horzcat(Channel_x, Channel_y);
% create output file
fSTAT = fullfile(rootP, 'STATtop.mat'); 


% read in the preprocessed data binary file
fid         = fopen(fsfile, 'rb'); % open for reading preprocessed data
if fid<3
    error('Could not open %s for reading.',fsfile);
end
s = dir(fsfile);
binary_size = s.bytes; %get the size of input file 


% read the spike times data
fprintf('Loading spike times... \n');
ss = h5read('D:\allKsOutput\Science\science.h5','/Dataset3/SpikeSample');


% read clustering data
fprintf('Loading clustering data... \n');
cluster = h5read('D:\allKsOutput\Science\science.h5','/Dataset3/SpikeCluster');
clusterChannels = h5read('D:\allKsOutput\Science\science.h5','/Dataset3/ClusterChannels'); %32 top channels used in each template
timeToCluster = horzcat(ss, cluster); %concatenate spike time, spike sample, cluster of the spike
sortedtimeToCluster = sortrows(timeToCluster,2); %sort according to cluster 


% set parameters
numTempSamp = 82; %number of samples extracted
nt0min = 20; %spike is aligned from 20 samples before the peak
nt0 = 61; %the total sample is 61



% Compute spike triggered average for templates 
labels = unique(cluster); %number of clusters
numCluster = numel(labels);
fprintf('Number of clusters: %d\n', numCluster);
[numChannel,~] = size(channelPosition);
% check if your snippet runs over the end of the file
maxSample = binary_size/(numChannel*2);
fprintf('Computing spike triggered average for templates... \n');
STAT = zeros(numCluster, numTempSamp, numChannel);
% select the top 32 channels as templates
[numtop, numall] = size(clusterChannels);
STATtop = zeros(numall, numTempSamp, numChannel);
for icluster = 1:numCluster
    frewind(fid);
    currentLabel = labels(icluster); %cluster number is not continuous
    spikes = sortedtimeToCluster(sortedtimeToCluster(:,2)==(currentLabel),:); %find the spikes in current cluster
    %allsnippet = zeros(length(spikes), numTempSamp, numChannel); %collection of snippet in one cluster
    summation = zeros(numTempSamp, numChannel, 'double');
    numSnippet = 0; %count of snippets
    fprintf('cluster %d, number of spike is %d\n', currentLabel, length(spikes));
    for jspike = 1:length(spikes)
        peak = spikes(jspike,1); %the peak sample number of a spike
        start = (peak-nt0min-(numTempSamp-nt0)-1)*numChannel*2; %get 41 samples before peak
        if start + numTempSamp < maxSample %make sure the snippet within bound
            fseek(fid, start, 'bof'); %get the waveform data at all channels of this spike 
            snippet = fread(fid, [numChannel numTempSamp], '*int16'); %read and reshape the datat to number of channels * number of samples in template
        %allsnippet(j,:,:) = snippet';
            snippet = snippet';
            summation = summation + double(snippet);
            numSnippet = numSnippet + 1;
        end
        if jspike == length(spikes)
            STAT(icluster,:,:) = summation/length(spikes);
        end
    end
    topChannel = sort(clusterChannels(:,currentLabel+1)); %get the top 32 channel numbers of a cluster 
    STATtop(currentLabel+1,(numTempSamp-nt0):end,topChannel) = STAT(icluster,(numTempSamp-nt0):end,topChannel);
end


save(fSTAT, 'STATtop'); %save stat as matlab matrix
fclose(fid);

fprintf('Complete! \n');



% Plot and compare individual template
KSTemps = h5read('D:\allKsOutput\Science\science.h5','/Dataset3/Templates');
testKS = KSTemps(1,:,:);
test2 = reshape(testKS,[],95);
figure(2)
plot(test2)

testSTAT = STATtop(1,:,:);
test1 = reshape(testSTAT,[],95);
figure(1)
plot(test1)

% testall = STAT(1,:,:);
% test3 = reshape(testall,[],95);
% figure(3)
% plot(test3)