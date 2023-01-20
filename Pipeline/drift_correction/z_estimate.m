function [diffZ,edges] = z_estimate(inputDir,fileName)
% This function estimated z_correction amount by the mode of histograms


% inputDir=input_struct.EMD_input_dir;
% fileName=matched_filename_corrected;

load(fullfile(inputDir,fileName));
hasMatch = ~isnan(f2_emd_ind);
f1_matched = f1(hasMatch,:);
f2_matched = f2(f2_emd_ind(hasMatch),:);
diffZ = f2_matched(:,2) - f1_matched(:,2);
binsize = 4;
edges = (-100:binsize:100);



