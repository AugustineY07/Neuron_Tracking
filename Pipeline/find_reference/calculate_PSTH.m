function calculate_PSTH(input_struct)
% This function computes rsp and PSTH from KS outputs with clusters selected based on 'good' choice, by default, all clusters are included
% Note all clusters, channels used here are converted to 1-based index

day = input_struct.day;
shank = input_struct.shank;
input_path = input_struct.rf_path;
output_path = input_struct.rf_path;
rootpath = input_struct.sorting_path;
input_name = input_struct.stimulus_data;
output_name = input_struct.psth_data;


stimulus = load(fullfile(input_path, input_name)).sorted_response; %stimulus_id:presention_time_s:presention_time_ms
day_label = load(fullfile(input_path, input_name)).day_label;
num_stim = 112;
nbins = 2002;


d = dir(rootpath); %list all files in directory
d = d([d.isdir]); %list sub-directories only
d = d(~ismember({d.name},{'.', '..'})); %exclude . and .. directories
smooth = gausswin(150, 3.5); %gaussian window used for PSTH smooth
for id = 1:day
    date = day_label(id);
    idx = find(contains({d.name},date));
    fday = fullfile(rootpath, d(idx).name);
    d1 = dir(fday);
    d1 = d1([d1.isdir]);
    d1 = d1(~ismember({d1.name},{'.', '..'}));
    for id1 = 1:shank
        fname = fullfile(fday, d1(id1).name, ['imec',num2str(id1-1),'_ks25']);
        csv = fullfile(fname, 'metrics.csv'); %summary data
        data = readtable(csv);
        fprintf('Finish data reading day %d shank %d. \n', id, id1);

        idx_good = [1:1:size(data,1)]; %include all clusters 

        good_clu(id,id1).clu = data.cluster_id(idx_good)+1; %zero-index to 1-index
        good_clu(id,id1).amp = data.amplitude(idx_good);
        good_clu(id,id1).ch = data.peak_channel(idx_good)+1; %zero-index to 1-index
        channel_position = readNPY(fullfile(fname, 'channel_positions.npy'));
        good_clu(id,id1).x = channel_position(good_clu(id,id1).ch,1);
        good_clu(id,id1).y = channel_position(good_clu(id,id1).ch,2);
        good_clu(id,id1).ISI = data.isi_viol(idx_good);
        fprintf('Good clusters found. \n');

        % compute rsp1, rsp2 - mean spike count stimulus to each of the num_stim natural images across trials (within 1 second after stimulus)
        spikeCluster = readNPY(fullfile(fname,'spike_clusters.npy'))+1;
        spikeSample = readNPY(fullfile(fname, 'spike_times.npy'));
        st = double(spikeSample)/30000;

        rspCount = zeros(length(good_clu(id,id1).ch),num_stim);
        for igroup = 1:size(stimulus,2) %presentation groups
            spikeByCluster = zeros(length(good_clu(id,id1).ch),num_stim); %hold the count of stimuli-triggered spikes in each cluster, dimension = good cluster*all stimulus
            if isempty(stimulus{id,igroup})
                continue
            else
                for istim = 1:num_stim
                    onset = stimulus{id,igroup}(istim,2); %time of stimulus presentation
                    offset = onset + 1; %1 second after stimulus presentation
                    start = find(st >= onset, 1); %get the index of first time point after presentation
                    ending = find(st <= offset, 1, 'last'); %get the index of last time point 1second after presentation
                    numSpike = ending - start;
                    all_sp_clu = spikeCluster(start:ending-1);
                    count = tabulate(all_sp_clu); %frequency of each value
                    if max(count(:,1)) <= max(good_clu(id,id1).clu)
                        count((max(count(:,1))+1):max(good_clu(id,id1).clu), :) = 0;
                    end
                    clu_count = count(good_clu(id,id1).clu,:); %only keep the valid clusters
                    spikeByCluster(:,istim) = clu_count(:,2);

                    %                     for ispike = 1:numSpike
                    %                         cluster = spikeCluster(start + ispike - 1); %associated cluster of the spike time above
                    %                         % find spike times in stimulus presentation window and add them up, added by cluster, not by channel
                    %                         if ismember(cluster,good_clu(id,id1).clu) == 1 %only include good cluster
                    %                             indChannel = find(good_clu(id,id1).clu == cluster); %find indexes of the cluster
                    %                             spikeByCluster(indChannel,istimulus) = spikeByCluster(indChannel,istimulus) + 1;
                    %                         end
                    %                     end
                end
                rspCount = rspCount + spikeByCluster; %holds stimulus sorted by index, every column is a stimuli
            end
            rspByTrial{id,id1,igroup} = spikeByCluster;
        end
        aveCount = rspCount/size(stimulus,2); %normalize rsp by number of trials
        rsp{id,id1} = aveCount;
        fprintf('Finish rsp calculation. \n');



        % compute rspPSTH1, rspPSTH2 - mean PSTH of each cluster across trials and stimulus (The PSTH is calculated with bins of 1ms)
        stms = st*1000; %spike time in ms
        PSTHByCluster = zeros(length(good_clu(id,id1).clu),2002); %from -0.5s to 1.5s on presentation
        for igroup = 1:size(stimulus,2)
            if isempty(stimulus{id,igroup})
                continue
            else
                for istimulus = 1:num_stim
                    tStimuli = stimulus{id,igroup}(istimulus,3); %stimulus presentation time
                    PSTHon = tStimuli - 501; %-0.5s
                    PSTHoff = tStimuli + 1501; %1.5s
                    %                     PSTHstart = find(stms >= PSTHon, 1); %find onset position in recording
                    %                     PSTHend = find(stms <= PSTHoff, 1, 'last');
                    %                     lenPSTH = PSTHend - PSTHstart; %number of spikes

                    start = find(stms >= PSTHon, 1); %get the index of first time point after presentation
                    ending = find(stms <= PSTHoff, 1, 'last');
                    all_sp_times = stms(start:ending-1); %all spike times in one stimulus
                    all_sp_cluster = spikeCluster(start:ending-1);

                    % count spikes by bins
                    bin_edges = [PSTHon:1:PSTHoff];
                    clu_edges = [good_clu(id,id1).clu', max(good_clu(id,id1).clu)+1]-0.1;
                    bin_clu_count = histcounts2(all_sp_cluster, all_sp_times, clu_edges, bin_edges);
                    PSTHByCluster = PSTHByCluster + bin_clu_count;
                end
            end
        end
        PSTH_mean{id,id1} = sum(PSTHByCluster,2)/(size(stimulus,2)*1.8*num_stim); %mean PSTH normalize by trial and time
        PSTHave = PSTHByCluster/(size(stimulus,2)*num_stim); %normalize by trial and image
        %PSTHcountByStimuli{id,id1} = PSTHave;
        PSTHsmooth = filter(smooth,1,PSTHave,[],2); %gaussian smooth along bins/columns
        PSTHsmoothed{id,id1} = PSTHsmooth(:,101:(end-100)); %remove 1 second at start and end
        fprintf('Finish PSTH calculation. \n');
    end
end

% pcount = 0;
% for iday = 1:size(PSTHsmoothed,1)
%     for ish = 1:size(PSTHsmoothed,2)
%         pcount = pcount + 1;
%         subplot(size(PSTHsmoothed,1),size(PSTHsmoothed,2),pcount)
%         plot(PSTHsmoothed{iday,ish}')
%     end
% end
% xlabel('Time(ms)')
% ylabel('Normalized Spike Count')
% [ax,h1]=suplabel('Days');
% set(h1,'FontSize',20)
% [ax,h2]=suplabel('Shanks','y');
% set(h2,'FontSize',20)

save(fullfile(output_path,output_name),'good_clu', 'rsp', 'rspByTrial', 'PSTHsmoothed', 'PSTH_mean');
