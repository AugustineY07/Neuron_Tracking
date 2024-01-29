function output = create_EMD_input(input, output, wf_metrics1, wf_metrics2, mwf1, mwf2, stage)

input_path = input.input_path;
EMD_input_dir = input.EMD_path;
data_path1 = input.data_path1;
data_path2 = input.data_path2;
KSLabel_name = input.KSLabel_name; % can be any tab delimited file of cluster# and labels. If bUseKSLabel = 1, only 'good' units will be included
v = input.validation;
all_z_mode = output.z_mode;

% output path
switch stage
    case 'pre'
        input_name = input.input_name;
    case 'post'
        input_name = input.input_name_post;
end
out_fullpath = fullfile(EMD_input_dir, input_name);


% 1 to include only units that pass KSlabel = good; keep set to 1 for now
bUseKSlabel = 1;

% get ks calls for each day

% load data
chan_pos = input.chan_pos;
mw1 = mwf1;
mw2 = mwf2;
wm1 = wf_metrics1;
switch stage
    case 'pre'
        wm2 = wf_metrics2;
    case 'post'
        wf_metrics2(:,10) = wf_metrics2(:,10) - all_z_mode; %z - mode
        wm2 = wf_metrics2;
end


% indices for:
% centX, centZ, fitY, pp, duration, fwhm, peak/trough ratio, pre-peak, z_footprint(z_sp)
metrics_ind = [9,10,11,1,3,4,5,6,8];


% Output: matlab file to read in EMD_main with pars.dataType = 'exp'
%  f1 -- nUnit x nDim for day 1 (e.g. x, z, pp_amp, pre_peak, duration, footprint)
%  f2 -- nUnit x nDim for day 2
%  If has validation: f2_same_ind -- index on day 1 of matches to day 2, Nan if no known match
%  metric_weights -- [nDim x 1] array of estimated error for each dimension
%  dim_mask = number of dimensions to include, e.g. just position or position + metrics

[npts_f1,~] = size(wm1);
[npts_f2,~] = size(wm2);
nDim = 9;
f1 = zeros([npts_f1,nDim]);
f2 = zeros([npts_f2,nDim]);


%load wf-fit locations
f1(:,1:3) = wm1(:,metrics_ind(1:3));
f2(:,1:3) = wm2(:,metrics_ind(1:3));
f1(:,4:9) = wm1(:,metrics_ind(4:9));
f2(:,4:9) = wm2(:,metrics_ind(4:9));


% check for units with no spikes (all mw entries = 0)
sum_unit_mw1 = sum(sum(mw1,[3]),[2]);
allzero_mw1 = (sum_unit_mw1 == 0);
sum_unit_mw2 = sum(sum(mw2,[3]),[2]);
allzero_mw2 = (sum_unit_mw2 == 0);

% Limit to particular shank
if isfield(input,'shank') && input.shank >= 0
    shank_center = 32 + input.shank*250;  % for 2.0 probes, chan_pos including offset
    onshank_mw1 = abs(wm1(:,metrics_ind(1))-shank_center) < 100;
    onshank_mw2 = abs(wm2(:,metrics_ind(1))-shank_center) < 100;
else
    onshank_mw1 = ones([npts_f1,1]);
    onshank_mw2 = ones([npts_f2,1]);
end

% get KSgood unit
if bUseKSlabel
    kscall1 = readKS2label(fullfile(input_path,data_path1,KSLabel_name),npts_f1);
    kscall2 = readKS2label(fullfile(input_path,data_path2,KSLabel_name),npts_f2);
    good_unit_1 = kscall1 & ~allzero_mw1 & onshank_mw1;
    good_unit_2 = kscall2 & ~allzero_mw2 & onshank_mw2;
else
    good_unit_1 = ~allzero_mw1 & onshank_mw1;
    good_unit_2 = ~allzero_mw2 & onshank_mw2;
end

fprintf('Number of units to include in match: %d, %d\n', sum(good_unit_1), sum(good_unit_2));

f1_good_orig = find(good_unit_1); % original 1-based unit labels for KSgood units
f2_good_orig = find(good_unit_2);
f1 = f1(good_unit_1,:);
f2 = f2(good_unit_2,:);
mw1 = mw1(good_unit_1,:,:);
mw2 = mw2(good_unit_2,:,:);
f1_labels = find(good_unit_1);
f2_labels = find(good_unit_2);


% If there is validation
if v == 1
    validation_path = input.validation_path;
    pairs = load(validation_path).ref;
    [npair, ~] = size(pairs);
    f2_same_ind = nan([npts_f1,1]);
    for i = 1:npair
        f2_same_ind(pairs(i,6)) = pairs(i,5);   % there's a row for each label, need to add 1 for matlab
    end
    f2_same_ind = nan([length(f1),1]);

    for i = 1:npair
        % check that the units in day 1 and day 2 are called as good
        f1_label = pairs(i,6);
        f2_label = pairs(i,5);

        if (ismember(f1_label,f1_good_orig) && ismember(f2_label, f2_good_orig))
            f1_new = find(f1_good_orig == f1_label);
            f2_new = find(f2_good_orig == f2_label);
            f2_same_ind(f1_new) = f2_new;
        end
    end
end




% will be stored
if ~exist(EMD_input_dir, 'dir')
    mkdir(EMD_input_dir);
end
if v == 1
    save(out_fullpath, 'f1', 'f2', 'f2_same_ind', 'mw1', 'mw2','chan_pos','f1_labels','f2_labels');
else
    save(out_fullpath, 'f1', 'f2','mw1', 'mw2','chan_pos','f1_labels','f2_labels');
end


end
