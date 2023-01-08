function [output_struct,pair_output] = find_reference(input_struct,stage,output_struct,pair_output,id,ish)

psth_data = input_struct.psth_data;
clu_data = input_struct.clu_data;
n = input_struct.n;
max_dist = input_struct.max_dist;
day = input_struct.day;
shank = input_struct.shank;
z_mode = pair_output.z_mode;
ref_path = input_struct.ref_path;


% load data
load(psth_data);
load(clu_data);
ref_pair_mode = [];
file_name = ['d1',num2str(id+1),'sh',num2str(ish),'_',stage,'_reference_4day_d1based.mat'];


numClu = simScore_sig{id,ish}; %simScore: day1*day2
fprintf('Looking for match between day %d and %d shank %d \n', id+1, 1, ish)
fprintf('There are %d and %d clusters \n', size(numClu,1), size(numClu,2))
dist = []; %distance to clu iclu

switch stage
    case 'pre'
        good_clu_KW;
    case 'post'
        good_clu_KW{id+1,ish}(:,4) = good_clu_KW{id+1,ish}(:,4) - z_mode; %drift correction
end

for iclu = 1:size(numClu,2) %num in day 2/col
    % rank simScore: pick n = neighbor tops; pick n = 8; set simScore > threshold
    [simRank,idxRank] = sort(numClu(:,iclu),'descend'); %rank idx in day 1/row
    top{id,ish,iclu}(:,1) = simRank(1:n,1); %highest n clusters
    top{id,ish,iclu}(:,2) = idxRank(1:n,1); %row number in the original(not sim_score sorted) order

    self = good_clu_KW{id+1,ish}(iclu,3:4); %location of current channel
    dist = sqrt((self(1)-good_clu_KW{1,ish}(:,3)).^2 + (self(2)-good_clu_KW{1,ish}(:,4)).^2);
    neighbor = find(dist <= max_dist); %all sorted index of 'neighbors'

    % check if neighbors contains at least 1 high simScore cluster and declare match:
    if any(ismember(neighbor,top{id,ish,iclu}(:,2)) == 1) %compare neighbor row idx with top idx in original order
        whichMatch = intersect(neighbor, top{id,ish,iclu}(:,2)); %find which one match, can be more than one
        fprintf('overlap rows are %d \n', whichMatch)

        switch length(whichMatch)
            case 1
                match{id,ish}(iclu,2) = whichMatch; %cluster row idx in good_clu_KW day 1 matches day 2
                match{id,ish}(iclu,4) = find(whichMatch==top{id,ish,iclu}(:,2)); %match rank in all simScore
                match{id,ish}(iclu,3) = top{id,ish,iclu}(match{id,ish}(iclu,4),1); %match simScore
            otherwise
                h_idx = arrayfun(@(x) find(top{id,ish,iclu}(:,2)==x,1),whichMatch); %idx in descending order
                match{id,ish}(iclu,3) = top{id,ish,iclu}(min(h_idx),1); %highest simScore
                match{id,ish}(iclu,2) = top{id,ish,iclu}(min(h_idx),2); %highest cluster row idx
                match{id,ish}(iclu,4) = min(h_idx); %highest rank
        end

        match{id,ish}(iclu,1) = iclu; %cluster row idx in day2
        match{id,ish}(iclu,5) = good_clu_KW{id+1,ish}(iclu,1); %cluster label in day2
        match{id,ish}(iclu,6) = good_clu_KW{1,ish}(match{id,ish}(iclu,2),1); %cluster label in day1
        match{id,ish}(iclu,7) = self(2) - good_clu_KW{1,ish}(match{id,ish}(iclu,2),4); % distance between the pair
    else
        match{id,ish}(iclu,:) = NaN;
    end
end


% exclude simscore < 1
match_large{id,ish} = match{id,ish}(match{id,ish}(:,3)>1,:);
match_large{id,ish} = sortrows(match_large{id,ish},3,'descend');

% find duplicates
[unique_clu,ua,~] = unique(match_large{id,ish}(:,6), 'stable'); %day1 clusters matched to
all_clu = match_large{id,ish}(:,6); %all day1 non-nan clusters match
duplicate = match_large{id,ish}(find(not(ismember(1:numel(all_clu),ua))),:); %all_clu idx that are duplicates
duplicate_temp = []; %get all duplicates
for idup = 1:size(duplicate,1)
    idx_dup = duplicate(idup,6); %duplicate day1 cluster
    dups = find(match_large{id,ish}(:,6) == idx_dup); %which rows contain duplicates
    duplicate_temp = [duplicate_temp; match_large{id,ish}(dups,:)];
end
duplicate_all{id,ish} = duplicate_temp;
if ~isempty(duplicate_all{id,ish})
    whichUnique = setdiff(match_large{id,ish}(:,1),duplicate_all{id,ish}(:,1)); %day2 idx that contains only unique matches
    nondup{id,ish} = match_large{id,ish}(find(ismember(match_large{id,ish}(:,1), whichUnique)),:);
else
    nondup{id,ish} = match_large{id,ish};
end

% ref_mode = 0;
% switch stage
%     case 'pre'
%         edges = (-100:binsize:100);
%         [num,~] = hist(nondup{id,ish}(:,7),edges);
%         values = max(num);
%         idx_mode = find(num==values); %all index
%         for iid = 1:length(idx_mode)
%             ref_mode = ref_mode + (edges(idx_mode(iid)-1)+edges(idx_mode(iid))) /(length(idx_mode)*2); %ave
%         end
%         ref_pair_mode(id+1,ish) = ref_mode;
% end
fprintf('Finish match! \n')


% generate output
ref_day1based = nondup{id,ish};
output_struct.KS = reshape(arrayfun(@(x) size(good_clu(x).clu,1), 1:numel(good_clu)), size(good_clu,1), size(good_clu,2)); 
output_struct.KW = cell2mat(cellfun(@length,good_clu_KW,'UniformOutput',false));


switch stage
    case 'pre'
        pair_output.match_pre = cell2mat(cellfun(@length,match,'UniformOutput',false)); %size(match,1);
        pair_output.ref_pre = cell2mat(cellfun(@length,nondup,'UniformOutput',false));
        pair_output.ref_filename_pre = file_name;
    case 'post'
        pair_output.match_post = cell2mat(cellfun(@length,match,'UniformOutput',false)); %size(match,1);
        pair_output.ref_post = cell2mat(cellfun(@length,nondup,'UniformOutput',false));
        pair_output.ref_filename_post = file_name;
end

% save data as files
if ~exist(ref_path, 'dir')
    mkdir(ref_path);
end
save(fullfile(ref_path, file_name), 'ref_day1based')
