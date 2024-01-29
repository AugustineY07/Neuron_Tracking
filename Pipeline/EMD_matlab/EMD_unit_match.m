% Match units using EMD
function output = EMD_unit_match(input,output,stage)


dim_mask = input.dim_mask;
dim_mask_physical = input.dim_mask_physical;
dim_mask_wf = input.dim_mask_wf;
l2_weight = input.l2_weights;
rootD = input.EMD_path;
v = input.validation;
xStep = input.xStep;
zStep = input.zStep;

if isfield(input,'diagDistCalc')
    diagDistCalc = input.diagDistCalc;
else
    diagDistCalc = false;
end

% output path
switch stage
    case 'pre'
        inputName = input.input_name;
        output_name = input.filename_pre; %matched file
    case 'post'
        inputName = input.input_name_post;
        output_name = input.filename_post; %matched file
end
out_fullpath = fullfile(rootD, output_name);

% load data: for experimental data, add mw1, mw2 chan_pos for these cluster
points = load(fullfile(rootD, inputName));
f1 = points.f1;
f2 = points.f2;
mw1 = points.mw1;
mw2 = points.mw2;
chan_pos = points.chan_pos;
f1_labels = points.f1_labels;
f2_labels = points.f2_labels;

% set weights = 1
w1 = ones(size(f1,1),1);%ones([length(f1),1])/length(f1);
w2 = ones(size(f2,1),1);%ones([length(f2),1])/length(f2);


% % test batch with #3 neurons 
% idx_f1 = [find(f1_labels==4),find(f1_labels==25),find(f1_labels==87)];
% idx_f2 = [find(f2_labels==4),find(f2_labels==120),find(f2_labels==151)];
% f1 = f1(idx_f1,:);
% f2 = f2(idx_f2,:);
% mw1 = mw1(idx_f1,:,:);
% mw2 = mw2(idx_f2,:,:);
% f1_labels = f1_labels(idx_f1,:);
% f2_labels = f2_labels(idx_f2,:);
% w1 = ones(size(f1,1),1);
% w2 = ones(size(f2,1),1);


% initialize distance mat
C = [];
Cp = [];
Cw = [];

switch stage
    case 'pre'
        [x, fval, C] = emd_nt(f1, f2, w1, w2, mw1, mw2, chan_pos, dim_mask, l2_weight, xStep, zStep, @weighted_gdf_nt);
        P = reshape(x,[size(f2,1),size(f1,1)]);
        if diagDistCalc
            [Cp,~] = gdm_nt(f1, f2, mw1, mw2, chan_pos, dim_mask_physical, l2_weight, xStep, zStep, @weighted_gdf_nt);
            [Cw,~] = gdm_nt(f1, f2, mw1, mw2, chan_pos, dim_mask_wf, l2_weight, xStep, zStep, @weighted_gdf_nt);
        else
            Cp = zeros(size(C));
            Cw = zeros(size(C));
        end
    case 'post'
        [x, fval, C] = emd_nt(f1, f2, w1, w2, mw1, mw2, chan_pos, dim_mask,l2_weight, xStep, zStep, @weighted_gdf_nt);
        P = reshape(x,[size(f2,1),size(f1,1)]);
        if diagDistCalc
            [Cp,~] = gdm_nt(f1, f2, mw1, mw2, chan_pos, dim_mask_physical, l2_weight, xStep, zStep, @weighted_gdf_nt);
            [Cw,~] = gdm_nt(f1, f2, mw1, mw2, chan_pos, dim_mask_wf, l2_weight, xStep, zStep, @weighted_gdf_nt);
        else
            Cp = zeros(size(C));
            Cw = zeros(size(C));
        end
end
cost = sum(C.*x);

C_unweighted_wf_dist = reshape(C,[size(f2,1),size(f1,1)]);
C_physical = reshape(Cp,[size(f2,1),size(f1,1)]);
C_wf = reshape(Cw,[size(f2,1),size(f1,1)]);


% summary matches or loop over the known matches (f2_same with values ~= NaN)
if v == 1
    nTP = 0;
    nFN = 0;
    nFP = 0;
    f2_same_ind = points.f2_same_ind;
    nPairs = sum(~isnan(f2_same_ind));
    pair_results = zeros([nPairs,9]);
    nf = 0; %counter for number of reference pairs encountered
end

np = numel(f1(:,1));
all_results = zeros([np,7]);
ip = 0; %counter for all pairs found

for nm = 1:np
    currCol = P(:,nm);
    [~, maxColInd] = max(currCol); % index of best match in f2

    % If unit nm in f1 was matched to a unit in f2, add it to the array of all_results
    % JIC note: moved calculation of physical and waveform distance to this
    % loop, so it is only performed for pairs.
    if max(currCol) > 0
        ip = ip + 1;
        all_results(ip,2) = f2_labels(maxColInd); %f2 label
        all_results(ip,3) = f1_labels(nm); %f1 label

        all_results(ip,4) = C_unweighted_wf_dist(maxColInd,nm); %loc+wf EMD distance

        [all_res_loc_dist, ~] = weighted_gdf_nt(f2(maxColInd, :), f1(nm, :), ...
            squeeze(mw2(maxColInd,:,:)), squeeze(mw1(nm,:,:)), chan_pos, ...
            dim_mask_physical, l2_weight, xStep, zStep);

        [all_res_wf_dist, ~] = weighted_gdf_nt(f2(maxColInd, :), f1(nm, :), ...
            squeeze(mw2(maxColInd,:,:)), squeeze(mw1(nm,:,:)), chan_pos, ...
            dim_mask_wf, l2_weight, xStep, zStep);

        all_results(ip,5) = all_res_loc_dist; % distance in physical space
        all_results(ip,6) = all_res_wf_dist; % waveform distance
        all_results(ip,7) = abs(f1(nm,2) - f2(maxColInd,2)); %z distance
        test_idx(nm) = maxColInd; %test if C is recorded correctly
    end  % if block for all pairs

    if v == 1
        if ~isnan(f2_same_ind(nm))
            % this f1 index (nm) had a reference match in f2
            nf = nf + 1;
            pair_results(nf,1) = f1_labels(nm);
            pair_results(nf,2) = f2_labels(f2_same_ind(nm)); % label of the ref match
            f2_ind = find(f2_labels==pair_results(nf,2)); % index of the ref match
            currRow = P(f2_ind,:);  % row that corresponds to the correct match for nm
            [~, maxRowInd] = max(currRow); % for a correct match, maxRowInd = nm

            if sum(currCol) == 0 && sum(currRow) == 0
                % no match called for either member this pair with any units
                nFN = nFN + 2;
                pair_results(nf,3) = -2;
                pair_results(nf,4) = -1;
                pair_results(nf,5) = -1;
                pair_results(nf,6) = 0;
                pair_results(nf,7) = 0;
                pair_results(nf,8) = 0;
                pair_results(nf,9) = 0;

            elseif sum(currCol) == 0 || sum(currRow) == 0
                % one member matched to a wrong unit,
                % the other unmatched
                nFP = nFP + 1;
                nFN = nFN + 1;
                pair_results(nf,3) = 1;
                if sum(currCol) == 0
                    % f1 has no pair, f2 = FP
                    pair_results(nf,4) = -1;
                    pair_results(nf,5) = f1_labels(maxRowInd);
                    % calculate distances for the incorrect paring of
                    % f2_ind to its EMD match.
                    pair_results(nf,6) = C_unweighted_wf_dist(f2_ind, maxRowInd); %loc+wf EMD distance
                    [loc_dist, ~] = weighted_gdf_nt(f2(f2_ind, :), f1(maxRowInd, :), ...
                        squeeze(mw2(f2_ind,:,:)), squeeze(mw1(maxRowInd,:,:)), chan_pos, ...
                        dim_mask_physical, l2_weight, xStep, zStep);
                    [wf_dist, ~] = weighted_gdf_nt(f2(f2_ind, :), f1(maxRowInd, :), ...
                        squeeze(mw2(f2_ind,:,:)), squeeze(mw1(maxRowInd,:,:)), chan_pos, ...
                        dim_mask_wf, l2_weight, xStep, zStep);
                    pair_results(nf,7) = loc_dist; % distance in physical space
                    pair_results(nf,8) = wf_dist; % waveform distance
                    pair_results(nf,9) = abs(f1(maxRowInd,2) - f2(f2_ind,2)); % z distance

%                 pair_results(nf,7) = C_physical(f2_ind,maxRowInd); %loc EMD distance
%                 pair_results(nf,8) = C_wf(f2_ind,maxRowInd); %wf EMD distance


                else
                    % f2 matched to this has no validated pair, so f1 = FP
                    % these are the distances already calculated for this
                    % pair in the all_results block above
                    pair_results(nf,4) = f2_labels(maxColInd);
                    pair_results(nf,5) = -1;
                    pair_results(nf,6) = C_unweighted_wf_dist(maxColInd,nm);
                    pair_results(nf,6) = all_res_loc_dist;
                    pair_results(nf,7) = all_res_wf_dist;
                    pair_results(nf,9) = abs(f1(nm,2) - f2(maxColInd,2));

%                     pair_results(nf,7) = C_physical(maxColInd,nm);
%                     pair_results(nf,8) = C_wf(maxColInd,nm); %nm,maxColInd
                   
                end
            elseif maxColInd == f2_same_ind(nm)
                % match to correct unit. distances are already calculated
                nTP = nTP + 2;
                pair_results(nf,4) = f2_labels(maxColInd);
                pair_results(nf,5) = f1_labels(maxRowInd);
                pair_results(nf,6) = C_unweighted_wf_dist(maxColInd,maxRowInd);
                pair_results(nf,7) = all_res_loc_dist;
                pair_results(nf,8) = all_res_wf_dist; 
                pair_results(nf,9) = abs(f1(maxRowInd,2) - f2(maxColInd,2));

            else
                % both units in a ref matched to other incorrect
                % units. Record the distances from te f2 ref unit to its incorrect EMD match 
                nFP = nFP + 2;
                pair_results(nf,3) = 2;
                pair_results(nf,4) = f2_labels(maxColInd);
                pair_results(nf,5) = f1_labels(maxRowInd);
                pair_results(nf,6) = C_unweighted_wf_dist(f2_ind, maxRowInd); %loc+wf EMD distance
                [loc_dist, ~] = weighted_gdf_nt(f2(f2_ind, :), f1(maxRowInd, :), ...
                    squeeze(mw2(f2_ind,:,:)), squeeze(mw1(maxRowInd,:,:)), chan_pos, ...
                    dim_mask_physical, l2_weight, xStep, zStep);
                [wf_dist, ~] = weighted_gdf_nt(f2(f2_ind, :), f1(maxRowInd, :), ...
                    squeeze(mw2(f2_ind,:,:)), squeeze(mw1(maxRowInd,:,:)), chan_pos, ...
                    dim_mask_wf, l2_weight, xStep, zStep);
                pair_results(nf,7) = loc_dist; % distance in physical space
                pair_results(nf,8) = wf_dist; % waveform distance
                pair_results(nf,9) = abs(f1(maxRowInd,2) - f2(f2_ind,2)); % z distance
% previous calculation, which gives distance from the incorrect match to f1
%  (maxColInd) to the incorrect match for f2(maxRowInd)
%                 pair_results(nf,6) = C_unweighted_wf_dist(maxColInd,maxRowInd);
%                 pair_results(nf,7) = C_physical(maxColInd,maxRowInd);
%                 pair_results(nf,8) = C_wf(maxColInd,maxRowInd); 
%                 pair_results(nf,9) = abs(f1(maxRowInd,2) - f2(maxColInd,2));
            end
        end
        
        % mark this row in all_results as having reference info
        pair_idx1 = find(pair_results(:,1) == all_results(ip,3));
        pair_idx2 = find(pair_results(:,2) == all_results(ip,2));
        if (pair_idx1 == pair_idx2 && pair_results(pair_idx1,3) == 0) %check if f1 and f2 label are in the same row
            all_results(ip,1) = 1; %has ref or not
        end
       
    end         % end of block for validation comparison


end
all_results = all_results(all_results(:,7)>0,:);


% count recovry rate if with validation
if v == 1
    if nPairs > 0
        fracCorr = nTP/(2*nPairs);
        fracFP = nFP/(2*nPairs);
        fracFN = nFN/(2*nPairs);
        fprintf('For experimental run: \n');
        fprintf('Frac correct %.3f\n', fracCorr);
        fprintf('Frac FP: %.3f\n', fracFP);
        fprintf('Frac FN: %.3f\n', fracFN);
    else
        fprintf('no pairs')
        fracCorr = [];
        fracFP = [];
        fracFN = [];
    end
end


% all units have pairs
f2_emd_ind = nan([length(f1),1]);

for nm = 1:length(f1)
    currCol = P(:,nm);
    if sum(currCol) ~= 0
        [~, maxColInd] = max(currCol);
        f2_emd_ind(nm) = maxColInd;
    end
end


% save results 
switch stage
    case 'pre'
        if v == 1 %with validation
            output.corr_pre = nTP/2;
            output.fp_pre = nFP/2;
            output.fn_pre = nFN/2;
            output.corrRate_pre = fracCorr;
            output.fpRate_pre = fracFP;
            output.fnRate_pre = fracFN;
            output.nPairs_pre = nPairs;
            output.pair_results_pre = pair_results;
        end
        output.C_pre = C;
        output.cost_pre = cost;  % note that this cost is the unweighted L2 * P
        output.P_pre = P;
        output.all_results_pre = all_results;
        output.x_pre = x;
    case 'post'
        if v == 1 %with validation
            output.corr_post = nTP/2;
            output.fp_post = nFP/2;
            output.fn_post = nFN/2;
            output.corrRate_post = fracCorr;
            output.fpRate_post = fracFP;
            output.fnRate_post = fracFN;
            output.nPairs_post = nPairs;
            output.pair_results_post = pair_results;
        end
        output.C_post = C;
        output.cost_post = cost;         % note that this cost is the unweighted L2 * P
        output.P_post = P;
        output.all_results_post = all_results;
        output.x_post = x;

        output.KSgood_f1 = size(f1,1);
        output.KSgood_f2 = size(f2,1);
        output.C_unweighted_wf_dist_post = C_unweighted_wf_dist;
        output.C_physical_post = C_physical;
        output.C_wf_post = C_wf;
        % output.L2 = C; removing, because this is a copy of C_unweighted_wf_dist
end

if v == 1
    save(out_fullpath,'f1','f2','f2_emd_ind','pair_results','f1_labels','f2_labels','all_results');
else
    save(out_fullpath,'f1','f2','f2_emd_ind','f1_labels','f2_labels','all_results');
end

end