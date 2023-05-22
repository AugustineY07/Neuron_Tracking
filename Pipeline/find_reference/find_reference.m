function [output_struct,pair_output] = find_reference(input_struct,stage,output_struct,pair_output,id,ish)

psth_data = input_struct.psth_data;
clu_data = input_struct.clu_data;
% n = input_struct.n;
max_dist = input_struct.max_dist;
day = input_struct.day;
shank = input_struct.shank;
z_mode = pair_output.z_mode;
ref_path = input_struct.rf_path;
out_path = input_struct.ref_path;
start_day = input_struct.start_day;

thre = max_dist;
% stage = 'pre';

% load data
load(fullfile(ref_path,psth_data));
load(fullfile(ref_path,clu_data));
% ref_pair_mode = [];

file_name = [input_struct.subject,'_d',num2str(start_day), 'd', num2str(id),'sh',num2str(ish),'_',stage,'_reference_4day','.mat'];
% switch compare
%     case '1-based'
%         file_name = [input_struct.subject,'d1',num2str(id+1),'sh',num2str(ish),'_',stage,'_reference_4day_d1based.mat'];
%     otherwise
%         file_name = [input_struct.subject,'d',num2str(id), 'd', num2str(id+1),'sh',num2str(ish),'_',stage,'_reference_4day_',compare,'.mat'];
% end

simScore = simScore_sig{id,ish}; %simScore: day1*day2
rspScore = sim_rsp_sig{id,ish}; 
psthScore = sim_PSTH_sig{id,ish}; 
fprintf('Looking for match between day %d and %d shank %d \n', id, start_day, ish)
fprintf('There are %d and %d clusters \n', size(simScore,1), size(simScore,2))
% dist = []; %distance to clu iclu

switch stage
    case 'pre'
        good_clu_KW;
    case 'post'
        good_clu_KW{id,ish}(:,4) = good_clu_KW{id,ish}(:,4) - z_mode; %drift correction of the later day
end





%% with physical constraints 
% for iclu = 1:size(simScore,2) %num in day 2/col
%     % rank simScore: pick n = neighbor tops; pick n = 8; set simScore > threshold
%     [simRank,idxRank] = sort(simScore(:,iclu),'descend'); %rank idx in day 1/row
%     top{id,ish,iclu}(:,1) = simRank(1:n,1); %highest n clusters
%     top{id,ish,iclu}(:,2) = idxRank(1:n,1); %row number in the original(not sim_score sorted) order
% 
%     self = good_clu_KW{id,ish}(iclu,3:4); %x,z location of current channel
%     dist = sqrt((self(1)-good_clu_KW{start_day,ish}(:,3)).^2 + (self(2)-good_clu_KW{start_day,ish}(:,4)).^2);
% %     switch compare
% %         case '1-based'
% %             dist = sqrt((self(1)-good_clu_KW{1,ish}(:,3)).^2 + (self(2)-good_clu_KW{1,ish}(:,4)).^2);
% %         case 'next'
% %             dist = sqrt((self(1)-good_clu_KW{id,ish}(:,3)).^2 + (self(2)-good_clu_KW{id,ish}(:,4)).^2);
% %     end
%     neighbor = find(dist <= max_dist); %all sorted index of 'neighbors'
% 
%     % check if neighbors contains at least 1 high simScore cluster and declare match:
%     if any(ismember(neighbor,top{id,ish,iclu}(:,2)) == 1) %compare neighbor row idx with top idx in original order
%         whichMatch = intersect(neighbor, top{id,ish,iclu}(:,2)); %find which one match, can be more than one
%         fprintf('overlap rows are %d \n', whichMatch)
% 
%         switch length(whichMatch)
%             case 1
%                 match{id,ish}(iclu,2) = whichMatch; %cluster row idx in good_clu_KW day 1 matches day 2
%                 match{id,ish}(iclu,4) = find(whichMatch==top{id,ish,iclu}(:,2)); %match rank in all simScore
%                 match{id,ish}(iclu,3) = top{id,ish,iclu}(match{id,ish}(iclu,4),1); %match simScore
%             otherwise
%                 h_idx = arrayfun(@(x) find(top{id,ish,iclu}(:,2)==x,1),whichMatch); %idx in descending order
%                 match{id,ish}(iclu,3) = top{id,ish,iclu}(min(h_idx),1); %highest simScore
%                 match{id,ish}(iclu,2) = top{id,ish,iclu}(min(h_idx),2); %highest cluster row idx
%                 match{id,ish}(iclu,4) = min(h_idx); %highest rank
%         end
% 
%         match{id,ish}(iclu,1) = iclu; %cluster row idx in day2
%         match{id,ish}(iclu,5) = good_clu_KW{id,ish}(iclu,1); %cluster label in day2
%         match{id,ish}(iclu,6) = good_clu_KW{start_day,ish}(match{id,ish}(iclu,2),1); %cluster label in day1
%         match{id,ish}(iclu,7) = self(2) - good_clu_KW{start_day,ish}(match{id,ish}(iclu,2),4); %distance between the pair
% %         switch compare
% %             case '1-based'
% %                 match{id,ish}(iclu,6) = good_clu_KW{1,ish}(match{id,ish}(iclu,2),1); 
% %                 match{id,ish}(iclu,7) = self(2) - good_clu_KW{1,ish}(match{id,ish}(iclu,2),4); 
% %             case 'next'
% %                 match{id,ish}(iclu,6) = good_clu_KW{id,ish}(match{id,ish}(iclu,2),1);
% %                 match{id,ish}(iclu,7) = self(2) - good_clu_KW{id,ish}(match{id,ish}(iclu,2),4); 
% %         end
%     else
%         match{id,ish}(iclu,:) = NaN;
%     end
% end
% 
% 
% % exclude simscore < 1
% match_large{id,ish} = match{id,ish}(match{id,ish}(:,3)>1,:);
% match_large{id,ish} = sortrows(match_large{id,ish},3,'descend');
% 
% % find duplicates
% [unique_clu,ua,~] = unique(match_large{id,ish}(:,6), 'stable'); %day1 clusters matched to
% all_clu = match_large{id,ish}(:,6); %all day1 non-nan clusters match
% duplicate = match_large{id,ish}(find(not(ismember(1:numel(all_clu),ua))),:); %all_clu idx that are duplicates
% duplicate_temp = []; %get all duplicates
% for idup = 1:size(duplicate,1)
%     idx_dup = duplicate(idup,6); %duplicate day1 cluster
%     dups = find(match_large{id,ish}(:,6) == idx_dup); %which rows contain duplicates
%     duplicate_temp = [duplicate_temp; match_large{id,ish}(dups,:)];
% end
% duplicate_all{id,ish} = duplicate_temp;
% if ~isempty(duplicate_all{id,ish})
%     whichUnique = setdiff(match_large{id,ish}(:,1),duplicate_all{id,ish}(:,1)); %day2 idx that contains only unique matches
%     nondup{id,ish} = match_large{id,ish}(find(ismember(match_large{id,ish}(:,1), whichUnique)),:);
% else
%     nondup{id,ish} = match_large{id,ish};
% end
% 
% % ref_mode = 0;
% % switch stage
% %     case 'pre'
% %         edges = (-100:binsize:100);
% %         [num,~] = hist(nondup{id,ish}(:,7),edges);
% %         values = max(num);
% %         idx_mode = find(num==values); %all index
% %         for iid = 1:length(idx_mode)
% %             ref_mode = ref_mode + (edges(idx_mode(iid)-1)+edges(idx_mode(iid))) /(length(idx_mode)*2); %ave
% %         end
% %         ref_pair_mode(id+1,ish) = ref_mode;
% % end
% fprintf('Finish match! \n')
% 
% 
% % generate output
% ref = nondup{id,ish};
% output_struct.KS = reshape(arrayfun(@(x) size(good_clu(x).clu,1), 1:numel(good_clu)), size(good_clu,1), size(good_clu,2)); 
% output_struct.KW = cell2mat(cellfun(@length,good_clu_KW,'UniformOutput',false));
% 
% 
% switch stage
%     case 'pre'
%         pair_output.match_pre = cell2mat(cellfun(@length,match,'UniformOutput',false)); %size(match,1);
%         pair_output.ref_pre = cell2mat(cellfun(@length,nondup,'UniformOutput',false));
%         pair_output.ref_filename_pre = file_name;
%     case 'post'
%         pair_output.match_post = cell2mat(cellfun(@length,match,'UniformOutput',false)); %size(match,1);
%         pair_output.ref_post = cell2mat(cellfun(@length,nondup,'UniformOutput',false));
%         pair_output.ref_filename_post = file_name;
% end
% 
% % save data as files
% if ~exist(ref_path, 'dir')
%     mkdir(ref_path);
% end
% save(fullfile(ref_path, file_name), 'ref')





%% without physical constraints 
assign_count = 1;
idx_next_round = [1:1:size(simScore,2)]; %initialize as all clusters
unit_simScore(idx_next_round,assign_count) = nan;
% stopping criteria: unassigned simscore <1, no unassigned unit
while ~isempty(idx_next_round) || ~all(unit_simScore(idx_next_round,assign_count-1)<1 & ~isnan(unit_simScore(idx_next_round,assign_count-1)))
    % for every round, reassign the next high simScore units to existing
    % match that: 1.>threshold; 2.has duplicate match
    unit_assignment(:,assign_count) = 0;
    unit_simScore(:,assign_count) = 0;
    for iclu = 1:length(idx_next_round) %num in day 2/col
        % rank simScore: set simScore > threshold
        [simRank,idxRank] = sort(simScore(:,idx_next_round(iclu)),'descend'); %rank idx in day 1/row
        rspRank = rspScore(idxRank,idx_next_round(iclu));
        psthRank = psthScore(idxRank,idx_next_round(iclu));

        % check simScore
        if simRank(assign_count) > 1
            % assign remaining highest simScore unit idx
            unit_assignment(idx_next_round(iclu),assign_count) = idxRank(assign_count);
            unit_simScore(idx_next_round(iclu),assign_count) = simRank(assign_count);
            unit_rank(idx_next_round(iclu),assign_count) = assign_count;
            unit_rsp(idx_next_round(iclu),assign_count) = rspRank(assign_count);
            unit_psth(idx_next_round(iclu),assign_count) = psthRank(assign_count);
        else
            % no match
            unit_assignment(idx_next_round(iclu),assign_count) = nan;
            unit_simScore(idx_next_round(iclu),assign_count) = nan;
            unit_rank(idx_next_round(iclu),assign_count) = nan;
            unit_rsp(idx_next_round(iclu),assign_count) = nan;
            unit_psth(idx_next_round(iclu),assign_count) = nan;
        end
    end

    % keep the assigned units
    if assign_count > 1
        unit_assignment(~idx_next_round_logical,assign_count) = unit_assignment(~idx_next_round_logical,assign_count-1);
        unit_simScore(~idx_next_round_logical,assign_count) = unit_simScore(~idx_next_round_logical,assign_count-1);
        unit_rank(~idx_next_round_logical,assign_count) = unit_rank(~idx_next_round_logical,assign_count-1);
        unit_rsp(~idx_next_round_logical,assign_count) = unit_rsp(~idx_next_round_logical,assign_count-1);
        unit_psth(~idx_next_round_logical,assign_count) = unit_psth(~idx_next_round_logical,assign_count-1);
    end

    % find duplicates: caveat -- unit assigned early but has a weak
    % simScore probably can't match to its next choice and ended up has no
    % match
    idx_nonan = find(~isnan(unit_assignment(:,assign_count)));
    exclude_nan = unit_assignment(idx_nonan,assign_count);
    [uniquevals, ~, ~] = unique(exclude_nan, 'stable');  % Stable keeps it in the same order
    [Ncount,edges] = histcounts(exclude_nan, [1:1:max(uniquevals)+1]);
    repeats = edges(Ncount > 1); %all duplicate units
    all_lost_idx = [];
    if ~isempty(repeats)
        for idup = 1:length(repeats)
            unit_repeat = repeats(idup);
            idx_repeat = find(unit_assignment(:,assign_count) == unit_repeat);
            [score,order] = sort(unit_simScore(idx_repeat),'descend');
            if score(1) > 1
                idx_lost_dup = idx_repeat(order(2:end)); %units with lower simscore, needs rematch
                all_lost_idx = [all_lost_idx idx_lost_dup'];
            end
        end
    end

    % check distance
    mask = find(isnan(unit_assignment(:,assign_count)));%idx of no-assignment
    unit_assignment(mask,assign_count) = 1; %assign an irrelevant unit
    dist(:,assign_count) = sqrt((good_clu_KW{id,ish}(:,3) - good_clu_KW{start_day,ish}(unit_assignment(:,assign_count),3)).^2 + (good_clu_KW{id,ish}(:,4) - good_clu_KW{start_day,ish}(unit_assignment(:,assign_count),4)).^2);
    dist(mask,assign_count) = nan; %get the nan back
    unit_assignment(mask,assign_count) = nan; 
    dist(all_lost_idx,assign_count) = 1000; %add duplicate unit back
    idx_next_round = find(dist(:,assign_count) > thre); %units need to rematch
    idx_next_round_logical = dist(:,assign_count) > thre;
    assign_count = assign_count + 1;
end 

    % declare match
    idx_nonan = find(~isnan(unit_assignment(:,assign_count-1)));
    match{id,ish}(:,1) = idx_nonan; %cluster row idx in day2
    match{id,ish}(:,2) = unit_assignment(idx_nonan,assign_count-1); 
    match{id,ish}(:,3) = unit_simScore(idx_nonan,assign_count-1); %highest simscore
    match{id,ish}(:,4) = unit_rank(idx_nonan,assign_count-1); %rank in all simscore
    match{id,ish}(:,5) = good_clu_KW{id,ish}(idx_nonan,1); %cluster label in day2
    match{id,ish}(:,6) = good_clu_KW{start_day,ish}(unit_assignment(idx_nonan,assign_count-1),1); %cluster label in day1
    match{id,ish}(:,7) = dist(idx_nonan,assign_count-1); %z distance between the pair
    match{id,ish}(:,8) = unit_rsp(idx_nonan,assign_count-1);
    match{id,ish}(:,9) = unit_psth(idx_nonan,assign_count-1);



% % exclude simscore < 1
% match_large{id,ish} = match{id,ish}(match{id,ish}(:,3)>1,:);
% match_large{id,ish} = sortrows(match_large{id,ish},3,'descend');
% 
% % find duplicates
% pair_unique = zeros(size(match_large{id,ish},1),size(match_large{id,ish},2));
% icount = 0;
% [unique_clu,ua,~] = unique(match_large{id,ish}(:,6), 'stable'); %day1 clusters
% for iu = 1:length(unique_clu)
%     numoccur = find(match_large{id,ish}(:,6) == unique_clu(iu));
%     if length(numoccur) > 1 %assign duplicated case with the largest simscore 
%         compare_sim = match_large{id,ish}(numoccur,3);
%         %compare_dist = abs(match_large{id,ish}(numoccur,7));
%         similar = max(compare_sim);
%         %short = min(compare_dist);
%         isim = find(compare_sim == similar);
%         %ishort = find(compare_dist == short);
%         if length(isim) == 1 
%             if abs(match_large{id,ish}(numoccur(isim),7)) <= thre
%                 icount = icount+1;
%                 pair_unique(icount,:) = match_large{id,ish}(numoccur(isim),:);
%             end
%         end
%     else
%         if abs(match_large{id,ish}(numoccur,7)) <= thre
%             icount = icount+1;
%             pair_unique(icount,:) = match_large{id,ish}(numoccur,:);
%         end
%     end
% end
% nondup{id,ish} = pair_unique(pair_unique(:,1)~=0,:); %drop 0 rows

fprintf('Finish match! \n')


% generate output
% ref = nondup{id,ish};
ref = match{id,ish};
output_struct.KS = reshape(arrayfun(@(x) size(good_clu(x).clu,1), 1:numel(good_clu)), size(good_clu,1), size(good_clu,2)); 
output_struct.KW = cell2mat(cellfun(@length,good_clu_KW,'UniformOutput',false));


switch stage
    case 'pre'
        pair_output.ref_pre = cell2mat(cellfun(@length,match,'UniformOutput',false));
        pair_output.ref_filename_pre = file_name;
    case 'post'
        pair_output.ref_post = cell2mat(cellfun(@length,match,'UniformOutput',false));
        pair_output.ref_filename_post = file_name;
end

% save data as files
if ~exist(out_path, 'dir')
    mkdir(out_path);
end
save(fullfile(out_path, file_name), 'ref')





















