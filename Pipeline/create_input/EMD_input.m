function EMD_input(input_struct, output_struct, pair_output, phydir1, phydir2, locdir1, locdir2, input_name,stage,start_day,id)

EMD_input_dir = input_struct.EMD_input_dir;
localization = input_struct.location;
ref_path = input_struct.ref_path;
rf_path = input_struct.rf_path;
subject = input_struct.subject;

% start_day = input_struct.start_day
switch stage
    case 'pre'
        fullpath = fullfile(ref_path, pair_output.ref_filename_pre);
        ref_name = pair_output.ref_filename_pre;
        all_z_mode = 0;
    case 'post'
        fullpath = fullfile(ref_path, pair_output.ref_filename_post);
        ref_name = pair_output.ref_filename_post;
        all_z_mode = pair_output.z_mode;
end

% 1 to include only units that pass KSlabel = good; keep set to 1 for now
bUseKSlabel = 1;

% get ks calls for each day
% fullpath = fullfile(phydir1,'cluster_KSLabel.tsv');
% kscall1 = ks_call;
% fullpath = fullfile(phydir2,'cluster_KSLabel.tsv');
% kscall2 = ks_call;
kscall1 = readKS2label(fullfile(phydir1,'cluster_KSLabel.tsv'));
kscall2 = readKS2label(fullfile(phydir2,'cluster_KSLabel.tsv'));

% get metrics 
chan_pos = readNPY(fullfile(phydir1,'channel_positions.npy'));
mw1 = readNPY(fullfile(phydir1, 'ksproc_mean_waveforms.npy'));
mw2 = readNPY(fullfile(phydir2, 'ksproc_mean_waveforms.npy'));

switch subject
    case 'AL031'
%         start_day = input_struct.start_day;
    load(fullfile(ref_path,ref_name));
    load(fullfile(rf_path,[subject,'_PSTH.mat']));
    %first day
    clu_label = ref(:,6); %get all clu label of ref units
    [~,same_idx1,~] = intersect(good_clu(start_day).clu, clu_label, 'stable');%get all clu idx in 'good_clu'
    y_loc = good_clu(start_day).y(same_idx1); %get y location of these channels
    y_lim = [min(y_loc) max(y_loc)];

    %second day
    clu_label2 = ref(:,5); %get all clu label of ref units
    [~,same_idx2,~] = intersect(good_clu(id).clu, clu_label2, 'stable');%get all clu idx in 'good_clu'
    y_loc2 = good_clu(id).y(same_idx2); %get y location of these channels
    y_lim2 = [min(y_loc2) max(y_loc2)];

end



wm1 = wave_metrics(mw1, chan_pos, 'pre',all_z_mode);
switch localization
    case 'Boussard'
        wm2 = wave_metrics(mw2, chan_pos, 'pre',all_z_mode);
    case 'wf'
        switch stage
            case 'pre'
                wm2 = wave_metrics(mw2, chan_pos, 'pre',all_z_mode);
            case 'post'
                z_mode = all_z_mode;
                wm2 = wave_metrics(mw2, chan_pos, 'post',all_z_mode);
        end
end

% indices for:
% centX, centZ, fitY, pp, duration, fwhm, peak/trough ratio, pre-peak, z_footprint(z_sp)
metrics_ind = [9,10,11,1,3,4,5,6,8];

% Output: matlab file to read in EMD_main with pars.dataType = 'exp'
%  f1 -- nUnit x nDim for day 1 (e.g. x, z, pp_amp, pre_peak, duration, footprint)
%  f2 -- nUnit x nDim for day 2
%  f2_same_ind -- index on day 1 of matches to day 2, Nan if no known match
%  metric_weights -- [nDim x 1] array of estimated error for each dimension
%  dim_mask = number of dimensions to include, e.g. just position or position + metrics
if ~isempty(fullpath)
    pairs = load(fullpath).ref; 
    [npair, ~] = size(pairs);
else
    npair = 0;
end


[npts_f1,~] = size(wm1);
[npts_f2,~] = size(wm2);
nDim = 9;
f1 = zeros([npts_f1,nDim]);
f2 = zeros([npts_f2,nDim]);

% right now loading all dimensions from the output of wave_metrics
% f1(:,1:2) = wm1(:,metrics_ind(1:2));
% f2(:,1:2) = wm2(:,metrics_ind(1:2));

switch localization
    case 'Boussard'
        % get centroids
        switch stage
            case 'pre'
                % get centroids
                cent1 = readNPY(fullfile(locdir1,'centXZY.npy'));
                cent2 = readNPY(fullfile(locdir2,'centXZY.npy'));
            case 'post'
                if exist([locdir1,'\centXZY_corrected.npy'],'file')
                    cent1 = readNPY(fullfile(locdir1,'centXZY_corrected.npy'));
                else
                    cent1 = readNPY(fullfile(locdir1,'centXZY.npy'));
                end
                cent2 = readNPY(fullfile(locdir2,'centXZY_corrected.npy'));
        end
        % load centroids from single spike localization
        f1(:,1:3) = cent1(:,:); % from single spike localizations
        f2(:,1:3) = cent2(:,:); % from single spike localizations

    case 'wf'
        switch subject
            case 'AL031'
                f1_idx = find(wm1(:,metrics_ind(2))<=y_lim(2) & wm1(:,metrics_ind(2))>=y_lim(1));
                f2_idx = find(wm2(:,metrics_ind(2))<=y_lim(2) & wm2(:,metrics_ind(2))>=y_lim(1));
                wm1 = wm1(f1_idx,:);
                wm2 = wm2(f2_idx,:);
                kscall1 = kscall1(f1_idx,:);
                kscall2 = kscall2(f2_idx,:);
                [npts_f1,~] = size(wm1);
                [npts_f2,~] = size(wm2);
                f1 = zeros([npts_f1,nDim]);
                f2 = zeros([npts_f2,nDim]);
                mw1 = mw1(f1_idx,:,:);
                mw2 = mw2(f2_idx,:,:);
                clu1 = readtable(fullfile(phydir1,'metrics.csv')).cluster_id+1;
                clu2 = readtable(fullfile(phydir2,'metrics.csv')).cluster_id+1;
                f1_idx_conversion(:,1) = [1:1:length(f1)]; %new idx after unit filtering
                [~,f1_idx_conversion(:,2),~] = intersect(clu1,f1_idx,'stable'); %old idx correponds to cluster idx 
                f1_idx_conversion(:,3) = f1_idx; %cluster label
                f2_idx_conversion(:,1) = [1:1:length(f2)];
                f2_idx_conversion(:,2) = f2_idx;
                [~,f2_idx_conversion(:,2),~] = intersect(clu2,f2_idx,'stable'); %cluster label
                f2_idx_conversion(:,3) = f2_idx;
        end
        %load wf-fit locations
        f1(:,1:3) = wm1(:,metrics_ind(1:3));
        f2(:,1:3) = wm2(:,metrics_ind(1:3));
end

f1(:,4:9) = wm1(:,metrics_ind(4:9));
f2(:,4:9) = wm2(:,metrics_ind(4:9));

f2_same_ind = nan([npts_f1,1]);
for i = 1:npair
    f2_same_ind(pairs(i,6)) = pairs(i,5);   % there's a row for each label, need to add 1 for matlab
end

% check for units with no spikes (all mw entries = 0)
sum_unit_mw1 = sum(sum(mw1,[3]),[2]);
allzero_mw1 = (sum_unit_mw1 == 0);
sum_unit_mw2 = sum(sum(mw2,[3]),[2]);
allzero_mw2 = (sum_unit_mw2 == 0);

if bUseKSlabel
    good_unit_1 = kscall1 & ~allzero_mw1;
    good_unit_2 = kscall2 & ~allzero_mw2;
else
    good_unit_1 = ~allzero_mw1;
    good_unit_2 = ~allzero_mw2;
end

f1_good_orig = find(good_unit_1); % original 1-based unit labels for good units
f2_good_orig = find(good_unit_2);
f1 = f1(good_unit_1,:);
f2 = f2(good_unit_2,:);
mw1 = mw1(good_unit_1,:,:);
mw2 = mw2(good_unit_2,:,:);
f2_same_ind = nan([length(f1),1]);
f1_labels = find(good_unit_1);
f2_labels = find(good_unit_2);

for i = 1:npair
    % check that the units in day 1 and day 2 are called as good
    f1_label = pairs(i,6);
    f2_label = pairs(i,5);
    switch subject
        case 'AL031'
            if (ismember(f1_label,f1_idx_conversion(f1_good_orig,3)) && ismember(f2_label, f2_idx_conversion(f2_good_orig,3)))
                f1_new = find(f1_idx_conversion(f1_good_orig,3) == f1_label);
                f2_new = find(f2_idx_conversion(f2_good_orig,3) == f2_label);
                f2_same_ind(f1_new) = f2_new;
            end
        otherwise
            if (ismember(f1_label,f1_good_orig) && ismember(f2_label, f2_good_orig))
                f1_new = find(f1_good_orig == f1_label);
                f2_new = find(f2_good_orig == f2_label);
                f2_same_ind(f1_new) = f2_new;
            end
    end
end



% will be stored in phydir 2
if ~exist(EMD_input_dir, 'dir')
    mkdir(EMD_input_dir);
end
switch subject
    case 'AL031'
        save(fullfile(EMD_input_dir,input_name), 'f1', 'f2', 'f2_same_ind', 'mw1', 'mw2','chan_pos','f1_labels','f2_labels','f1_idx_conversion','f2_idx_conversion');
    otherwise
        save(fullfile(EMD_input_dir,input_name), 'f1', 'f2', 'f2_same_ind', 'mw1', 'mw2','chan_pos','f1_labels','f2_labels');
end
end








%------------Helper function------------

function ks_call = readKS2label(fullpath)

fid = fopen(fullpath,'r');
nUnit = 0;
%allocate excess space to hold ks calls
ks_call = zeros([5000,1], 'logical');

% read line (this will be the header)
tline = fgetl(fid);
% read another line (first unit entry
tline = fgetl(fid);
while ischar(tline)
    nUnit = nUnit + 1;
    ln = split(tline);
    call_str = ln{2};
    if strcmp(call_str,'good')
        ks_call(nUnit) = 1;
    end
    tline = fgetl(fid);
end

ks_call = ks_call(1:nUnit);
fprintf('%d out of %d units called good\n', sum(ks_call), nUnit);
end

