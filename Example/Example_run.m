% An example of tracking units in animal AL032 shank 1 

clear;

%----------add packages----------
% EDIT THESE PATHS FOR YOUR SYSTEM
addpath(genpath('C:\Users\colonellj\Documents\matlab_ephys\AY_match_current\Neuron_Tracking')) %Add this repo to MATLAB path 
addpath(genpath('C:\Users\colonellj\Documents\npy-matlab-master')) %Add npy-matlab to MATLAB path

%----------define parameter and path----------
%parent directory of the data (see Example directory in this repo for organizaiont) -- EDIT FOR YOUR SYSTEM 
input.input_path = 'C:\Users\colonellj\Documents\matlab_ephys\AY_match_current\Neuron_Tracking\Example'; 
input.EMD_path = input.input_path;
% Input data
input.fs = 30000; %acquisition rate
input.ts = 82; %wf time samples
input.threshold = 10;    % distance threshold for calling matches as real
input.l2_weights = 1500; % weight for waveform difference, see comments in Pipeline\EMD_matlab\weighted_gdf_nt
input.dim_mask = logical([1,1,1,0,0,0,0,0,0,1]); % used to calculate full distance metric
input.dim_mask_physical = logical([1,1,1,0,0,0,0,0,0,0]); % used to calculate 'position only' part of distance metric
input.dim_mask_wf = logical([0,0,0,0,0,0,0,0,0,1]); % used to calculate 'waveform only' part of distance metric
input.chan_pos_name = 'channel_positions.npy';
input.wf_name = 'ksproc_mean_waveforms.npy';
input.KSLabel_name = 'cluster_KSLabel.tsv';
input.validation = 0;

numData = 5; 



%----------Unit tracking----------
% Find match of all datasets (default: day n and day n+1, can be changed to track between non-consecutive datasets)
for id = 1:numData-1
    input.data_path1 = ['D',num2str(id)]; % first dataset, CHANGE to track non-consecutively 
    input.data_path2 = ['D',num2str(id+1)]; % second dataset, CHANGE to track non-consecutively
    input.result_path = fullfile(input.input_path,['result',num2str(id),num2str(id+1)]); %result directory 
    input.input_name = ['input',num2str(id),'.mat']; 
    input.input_name_post = ['input_post',num2str(id),'.mat']; 
    input.filename_pre = ['EMD_pre',num2str(id),'.mat']; 
    input.filename_post = ['EMD_post',num2str(id),'.mat']; 
    chan_pos = readNPY(fullfile(input.input_path, input.data_path1, 'channel_positions.npy'));
    mwf1 = readNPY(fullfile(input.input_path, input.data_path1, 'ksproc_mean_waveforms.npy'));
    mwf2 = readNPY(fullfile(input.input_path, input.data_path2, 'ksproc_mean_waveforms.npy'));
    NT_main(input,chan_pos,mwf1,mwf2);
end

% Find chains
for id = 1:numData-1 % Load data
    all_input(id) = load(fullfile(input.input_path,['result',num2str(id),num2str(id+1)], "Input.mat")); 
    all_output(id) = load(fullfile(input.input_path,['result',num2str(id),num2str(id+1)],'Output.mat'));
end
[chain_all,z_loc,len] = chain_summary(all_input,all_output,numData,input.input_path);



%----------Plot chains of interest (waveform, firing rate, original location, drift-corrected location, L2)----------
full_chain = chain_all(len == numData,:); %find chains with length across all datasets
[L2_weight,fr_all,fr_change,x_loc_all,z_loc_all] = chain_stats(all_input,all_output,full_chain,numData,input.input_path);

numChain = size(full_chain,1);
ichain = 1; %which chain to plot, please enter a number between 1 and numChain as input, NEED CHANGE  

figure()
for id = 1:numData-1
    % plot waveform
    plot_wf(all_input, full_chain, L2_weight, chan_pos, numData, ichain, id);
end
% plot firing rate
plot_fr(fr_all, fr_change, numData, ichain);
% plot location
plot_loc(all_input,x_loc_all,z_loc_all, chan_pos, numData,ichain)


