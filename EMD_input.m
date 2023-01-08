function EMD_input(input_struct, output_struct, pair_output, phydir1, phydir2, locdir1, locdir2, input_name,stage,ish, id)

EMD_input_dir = input_struct.EMD_input_dir;
localization = input_struct.location;
all_z_mode = pair_output.z_mode;
switch stage
    case 'pre'
        fullpath = pair_output.ref_filename_pre;
    case 'post'
        fullpath = pair_output.ref_filename_post;
end

% 1 to include only units that pass KSlabel = good; keep set to 1 for now
bUseKSlabel = 1;

% get ks calls for each day
kscall1 = readKS2label(fullfile(phydir1,'cluster_KSLabel.tsv'));
kscall2 = readKS2label(fullfile(phydir2,'cluster_KSLabel.tsv'));

% get metrics
chan_pos = readNPY(fullfile(phydir1,'channel_positions.npy'));
mw1 = readNPY(fullfile(phydir1, 'ksproc_mean_waveforms.npy'));
mw2 = readNPY(fullfile(phydir2, 'ksproc_mean_waveforms.npy'));


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
                wm2 = wave_metrics(mw2, chan_pos, 'post',z_mode);
        end
end

% indices for:
% centX, centZ, Y, pp, duration, fwhm, peak/trough, pre-peak, z footprint
metrics_ind = [9,10,11,1,3,4,5,6,8];

% Output: matlab file to read in EMD_main with pars.dataType = 'exp'
%  f1 -- nUnit x nDim for day 1 (e.g. x, z, pp_amp, pre_peak, duration, footprint)
%  f2 -- nUnit x nDim for day 2
%  f2_same_ind -- index on day 1 of matches to day 2, Nan if no known match
%  metric_weights -- [nDim x 1] array of estimated error for each dimension
%  dim_mask = number of dimensions to include, e.g. just position or position + metrics
if ~isempty(fullpath)
    pairs = load(fullpath).ref_day1based;
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
    if (ismember(f1_label,f1_good_orig) && ismember(f2_label, f2_good_orig))
        f1_new = find(f1_good_orig == f1_label);
        f2_new = find(f2_good_orig == f2_label);
        f2_same_ind(f1_new) = f2_new;
    end
end

% will be stored in phydir 2
if ~exist(EMD_input_dir, 'dir')
    mkdir(EMD_input_dir);
end
save(fullfile(EMD_input_dir,input_name), 'f1', 'f2', 'f2_same_ind', 'mw1', 'mw2','chan_pos','f1_labels','f2_labels');
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

