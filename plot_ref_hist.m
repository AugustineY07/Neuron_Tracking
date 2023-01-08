function plot_ref_hist(input_struct, pair_output, stage, id, ish)

opt = input_struct.plot_mask(1);
ref_path = input_struct.ref_path;
switch stage
    case 'pre'
        file_name = pair_output.ref_filename_pre;
    case 'post'
        file_name = pair_output.ref_filename_post;
end

load(fullfile(ref_path, file_name))


switch opt
    case 1
        % Plot ref histograms
        binsize = 4;
        edges = (-100:binsize:100);
        figure()
        histogram(ref_day1based(:,7),edges);
        title(sprintf('%s-correction all reference distance histogram day %d-1 shank %d', stage, id+1, ish))
end
