% Match units using EMD
function output = EMD_unit_match(input,output,stage)

dim_mask = input.dim_mask;
dim_mask_physical = input.dim_mask_physical;
dim_mask_wf = input.dim_mask_wf;
l2_weight = input.l2_weights;
rootD = input.EMD_path;
v = input.validation;

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

% initialize distance mat
C = [];
Cp = [];
Cw = [];

switch stage
    case 'pre'
        [x, fval, L2] = emd_nt(f1, f2, w1, w2, mw1, mw2, chan_pos, dim_mask, l2_weight, @weighted_gdf_nt);
        P = reshape(x,[length(f2),length(f1)]);
        C = gdm_nt(f1, f2, mw1, mw2, chan_pos, dim_mask, l2_weight, @weighted_gdf_nt); % Distance matrix
        Cp = gdm_nt(f1, f2, mw1, mw2, chan_pos, dim_mask_physical, l2_weight, @weighted_gdf_nt);
        Cw = gdm_nt(f1, f2, mw1, mw2, chan_pos, dim_mask_wf, l2_weight, @weighted_gdf_nt);
    case 'post'
        [x, fval, L2] = emd_nt(f1, f2, w1, w2, mw1, mw2, chan_pos, dim_mask,l2_weight, @weighted_gdf_nt);
        P = reshape(x,[length(f2),length(f1)]);
        C = gdm_nt(f1, f2, mw1, mw2, chan_pos, dim_mask, l2_weight, @weighted_gdf_nt);
        Cp = gdm_nt(f1, f2, mw1, mw2, chan_pos, dim_mask_physical, l2_weight, @weighted_gdf_nt);
        Cw = gdm_nt(f1, f2, mw1, mw2, chan_pos, dim_mask_wf, l2_weight, @weighted_gdf_nt);
end
cost = sum(C.*x);
C_dist = reshape(C,[length(f1),length(f2)]);
C_physical = reshape(Cp,[length(f1),length(f2)]);
C_wf = reshape(Cw,[length(f1),length(f2)]);




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
    [~, maxColInd] = max(currCol);

    if v == 1
        if ~isnan(f2_same_ind(nm))
            nf = nf + 1;
            pair_results(nf,1) = f1_labels(nm);
            pair_results(nf,2) = f2_labels(f2_same_ind(nm));
            f2_ind = find(f2_labels==pair_results(nf,2));
            currRow = P(f2_ind,:);
            [~, maxRowInd] = max(currRow);

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
                    pair_results(nf,6) = C_dist(maxRowInd,f2_ind); %loc+wf EMD distance
                    pair_results(nf,7) = C_physical(maxRowInd,f2_ind); %loc EMD distance
                    pair_results(nf,8) = C_wf(maxRowInd,f2_ind); %wf EMD distance
                    pair_results(nf,9) = abs(f1(maxRowInd,2) - f2(f2_ind,2)); %loc EMD z_distance

                else
                    % f2 has no pair, f1 = FP
                    pair_results(nf,4) = f2_labels(maxColInd);
                    pair_results(nf,5) = -1;
                    pair_results(nf,6) = C_dist(nm,maxColInd);
                    pair_results(nf,7) = C_physical(nm,maxColInd);
                    pair_results(nf,8) = C_wf(nm,maxColInd); %nm,maxColInd
                    pair_results(nf,9) = abs(f1(nm,2) - f2(maxColInd,2));

                end
            elseif maxColInd == f2_same_ind(nm)
                % match to correct unit
                nTP = nTP + 2;
                pair_results(nf,4) = f2_labels(maxColInd);
                pair_results(nf,5) = f1_labels(maxRowInd);
                pair_results(nf,6) = C_dist(maxRowInd,maxColInd);
                pair_results(nf,7) = C_physical(maxRowInd,maxColInd); %C_physical(maxRowInd,maxColInd)
                pair_results(nf,8) = C_wf(maxRowInd,maxColInd);
                pair_results(nf,9) = abs(f1(maxRowInd,2) - f2(maxColInd,2));

            else
                % both units matched to other incorrect
                % units
                nFP = nFP + 2;
                pair_results(nf,3) = 2;
                pair_results(nf,4) = f2_labels(maxColInd);
                pair_results(nf,5) = f1_labels(maxRowInd);
                pair_results(nf,6) = C_dist(maxRowInd,maxColInd);
                pair_results(nf,7) = C_physical(maxRowInd,maxColInd);
                pair_results(nf,8) = C_wf(maxRowInd,maxColInd);
                pair_results(nf,9) = abs(f1(maxRowInd,2) - f2(maxColInd,2));
            end
        end
    end

    % all results
    if max(currCol) > 0
        ip = ip + 1;
        all_results(ip,2) = f2_labels(maxColInd); %f2 label
        all_results(ip,3) = f1_labels(nm); %f1 label
        if v == 1
            pair_idx1 = find(pair_results(:,1) == all_results(ip,3));
            pair_idx2 = find(pair_results(:,2) == all_results(ip,2));
            if (pair_idx1 == pair_idx2 & pair_results(pair_idx1,3) == 0) %check if f1 and f2 label are in the same row
                all_results(ip,1) = 1; %has ref or not
            end
        end
        all_results(ip,4) = C_dist(nm,maxColInd); %loc+wf EMD distance
        all_results(ip,5) = C_physical(nm,maxColInd); %loc EMD distance
        all_results(ip,6) = C_wf(nm,maxColInd); %wf EMD distance
        all_results(ip,7) = abs(f1(nm,2) - f2(maxColInd,2)); %z distance
    end
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
        output.cost_pre = cost;
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
        output.cost_post = cost;
        output.P_post = P;
        output.all_results_post = all_results;
        output.x_post = x;

        output.KSgood_f1 = size(f1,1);
        output.KSgood_f2 = size(f2,1);
        output.C_dist_post = C_dist;
        output.C_physical_post = C_physical;
        output.C_wf_post = C_wf;
        output.L2 = L2;
end

if v == 1
    save(out_fullpath,'f1','f2','f2_emd_ind','pair_results','f1_labels','f2_labels','all_results');
else
    save(out_fullpath,'f1','f2','f2_emd_ind','f1_labels','f2_labels','all_results');
end

end