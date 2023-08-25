clear all

fig_path = 'C:\Users\labadmin\Desktop\Neuron Tracking Pipeline\User version\Github Figure reproduce\F5&6'; %NEED CHANGE to local path
load(fullfile(fig_path,'F5_data.mat'))

subject = {'AL031','AL032','AL036'};%,'AL032','AL036'
for is = 1:length(subject)
    switch subject{is}
        case 'AL036'
            day = 6;
            shank = 4;
            dur = [1 2 4 10 16 24 36];
            c = [0.0078 0 1; 0 0.4609 1;0 0.9297 1];
        case 'AL031'
            day = 4;
            shank = 1;
            dur = [1 2 11 36 45];
            c = [0 0.5 0; 0.4667 0.7333 0; 0.8588 0.9294 0];
        case 'AL032'
            day = 4;
            shank = 4;
            dur = [1 2 13 23 48];
            c = [1 0.5391 0; 1 0.75 0;1 0.906 0];
    end

    %xticks
    for ir = 1:day+1
        xlb{ir} = ['day ',num2str(dur(ir))];
    end

    ref_dur = [];
    putative_dur = [];
    mixed_dur = [];
    for ish = 1:shank
        % read data
        chains = track_all{is,ish}; %all chains
        ref_condition = ref_idx{is,ish};
        chains_full = chain_full{is,ish}; %full chains
        ref_condition_full = with_ref_full{is,ish};

        % compute chain length
        for ichain = 1:size(chains,1)
            idx_start = find(chains(ichain,:) ~= 0,1,'first'); %first and last nonzero value
            idx_end = find(chains(ichain,:) ~= 0,1,'last');
            chain_dur{is,ish}(ichain,1) = dur(idx_end) - dur(idx_start);

            %assign ref_condition to all chains
            idx_condition_start = find(ref_condition(ichain,:) ~= 0,1,'first'); %first and last nonzero value in the conditions
            idx_condition_end = find(ref_condition(ichain,:) ~= 0,1,'last'); 
            thisCon = ref_condition(ichain,idx_condition_start:idx_condition_end);
            if all(thisCon(1,:) == 2)
                condition_idx{is,ish}(ichain,1) = 1;
            elseif all(thisCon(1,:) == 1)
                condition_idx{is,ish}(ichain,1) = 2;
            else
                condition_idx{is,ish}(ichain,1) = 3;
            end
        end %end of chain

        % split into three types and merge shanks
        ref_index = find(condition_idx{is,ish}==1);
        ref_dur = [ref_dur; chain_dur{is,ish}(ref_index,:)];
        putative_idx = find(condition_idx{is,ish}==2);
        putative_dur = [putative_dur; chain_dur{is,ish}(putative_idx,:)];
        mixed_idx = find(condition_idx{is,ish}==3);
        mixed_dur = [mixed_dur; chain_dur{is,ish}(mixed_idx,:)];
    end %end of shank

    % survival plot by condition
    [day_count_ref,day_edge_ref] = histcounts(ref_dur,'BinLimits',[0.5,47.5]); %count each duration, 1 increment
    nz_idx_ref = find(day_count_ref~=0);
    nz_day_ref = day_edge_ref(nz_idx_ref); %nonzero durations
    nz_count_ref = day_count_ref(nz_idx_ref); %nonzero counts
    for isum = 1:length(nz_idx_ref) %cumulative num of chains 
        survival_chain_count_ref(isum) = sum(nz_count_ref(1:length(nz_count_ref)-isum+1));
    end
    nz_day_ref(2:length(nz_day_ref)+1) = nz_day_ref+0.5; %shift right to connect the lines
    nz_day_ref(1) = 0;
    survival_chain_count_ref(2:length(survival_chain_count_ref)+1) = survival_chain_count_ref;
    stairs(nz_day_ref,survival_chain_count_ref,'LineWidth',3,'Marker','*','color',c(1,:)); hold on %ref

    [day_count_put,day_edge_put] = histcounts(putative_dur,'BinLimits',[0.5,47.5]); 
    nz_idx_put = find(day_count_put~=0);
    nz_day_put = day_edge_put(nz_idx_put); 
    nz_count_put = day_count_put(nz_idx_put); 
    for isum = 1:length(nz_idx_put) 
        survival_chain_count_put(isum) = sum(nz_count_put(1:length(nz_count_put)-isum+1));
    end
    nz_day_put(2:length(nz_day_put)+1) = nz_day_put+0.5; %shift right to connect the lines
    nz_day_put(1) = 0;
    survival_chain_count_put(2:length(survival_chain_count_put)+1) = survival_chain_count_put;
    stairs(nz_day_put,survival_chain_count_put,':','LineWidth',3,'Marker','*','color',c(2,:)); hold on %putative

    [day_count_mixed,day_edge_mixed] = histcounts(mixed_dur,'BinLimits',[0.5,47.5]); 
    nz_idx_mixed = find(day_count_mixed~=0);
    nz_day_mixed = day_edge_mixed(nz_idx_mixed); 
    nz_count_mixed = day_count_mixed(nz_idx_mixed); 
    for isum = 1:length(nz_idx_mixed) 
        survival_chain_count_mixed(isum) = sum(nz_count_mixed(1:length(nz_count_mixed)-isum+1));
    end
    nz_day_mixed(2:length(nz_day_mixed)+1) = nz_day_mixed+0.5; %shift right to connect the lines
    nz_day_mixed(1) = 0;
    survival_chain_count_mixed(2:length(survival_chain_count_mixed)+1) = survival_chain_count_mixed;
    stairs(nz_day_mixed,survival_chain_count_mixed,'--','LineWidth',3,'Marker','*','color',c(3,:)); hold on %mixed
end %end of subject


legend('AL031 ref', 'AL031 putative', 'AL031 mixed','AL032 ref', 'AL032 putative', 'AL032 mixed','AL036 ref', 'AL036 putative', 'AL036 mixed', 'Location', 'northeast','FontSize',18); hold on;
xlim([0 50])
ylim([0 80])
xlabel('Tracking duration (\Deltadays)','FontSize',18,'FontWeight','Bold','Color','k')
ylabel('Cumulative number of units tracked','FontSize',18,'FontWeight','Bold','Color','k')
title('Summary of duration of neuron tracked across all subjects','FontSize',18,'FontWeight','Bold','Color','k')

ax = gca;
ax.FontSize = 16; %tick font
ax.Box = 'off'; %remove tick box
set(ax,'TickDir','out'); %tickmark towards outside


