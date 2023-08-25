% Calculate recovery rate and accuracy if has validation
function [eval] = acc(varargin)

if isempty(varargin)
    [rezName, rezPath] = uigetfile('*.mat','Select validation table for the two datasets');
    truth = load(fullfile(rezPath, rezName)).label; %'C:\Users\labadmin\Desktop\Neuron Tracking Pipeline\User version\Supplementary\truth.mat'

    [mmName, mmPath] = uigetfile('*.mat','Select matching output for the two datasets(example: in result12)'); %'C:\Users\labadmin\Desktop\Neuron Tracking Pipeline\User version\result12\Output.mat'
    pair_output = load(fullfile(mmPath, mmName)).output;
end

all_match = pair_output.all_results_post;
threshold_match = pair_output.results_wth;

% recovery rate
[same_clu,same_idx1,same_idx2] = intersect(truth(:,2), all_match(:,3),'stable');
label_matched = truth(same_idx1,:);
all_matched = all_match(same_idx2,:);
all_matched = all_match(same_idx2,2:3);
[q, all_idx] = ismember(label_matched, all_matched, 'rows');
eval(1) = nnz(all_idx)/length(all_idx);

% accuracy
[same_clu_th,same_idx1_th,same_idx2_th] = intersect(truth(:,2), threshold_match(:,3),'stable');
label_matched_th = truth(same_idx1_th,:);
threshold_matched = threshold_match(same_idx2_th,2:3);
[~, th_idx] = ismember(label_matched_th, threshold_matched, 'rows'); 
eval(2) = nnz(th_idx)/length(th_idx);

end