% An example of tracking units in animal AL032 shank 1 

%----------add packages----------
addpath(genpath('C:\Users\labadmin\Desktop\Neuron Tracking Pipeline\User version')) %NEED CHANGE 
addpath(genpath('D:\Data\Pipeline\npy\npy-matlab')) %path to your npy directory

%----------define parameter and path----------
input.input_path = 'C:\Users\labadmin\Desktop\Neuron Tracking Pipeline\User version\'; %main directory, NEED CHANGE 
input.EMD_path = fullfile(input.input_path,'EMD_input\'); %EMD input directory, NEED CHANGE 
% Input data
input.fs = 30000; %acquisition rate, NEED to match your dataset!
input.ts = 82; %wf time samples
input.l2_weights = 1500;
input.threshold = 10; %z distance threshold for matched units, 10um recommended, can change if needed
input.validation = 0; %no reference data, set to 1 if you have ground truth
input.xStep = 32; %space between columns of sites, um (NP 2.0 = 32, can find in the channel map)
input.zStep = 15; %space between rows of sites, um (NP 2.0 = 15, can find in the channel map)
input.dim_mask = logical([1,1,1,0,0,0,0,0,0,1]); %default = x,z,y position, waveform distance
input.dim_mask_physical = logical([1,1,1,0,0,0,0,0,0,0]);
input.dim_mask_wf = logical([0,0,0,0,0,0,0,0,0,1]);
input.chan_pos_name = 'channel_positions.npy'; %NEED to match your channel map file name!
input.wf_name = 'ksproc_mean_waveforms.npy'; %NEED to match your waveform file name!
input.KSLabel_name = 'cluster_KSLabel.tsv'; %NEED to match your KSlabel file name!

numData = 5; %Number of datasets to match



%----------Unit tracking----------
if exist(input.EMD_path, 'dir') == 0
    mkdir(input.EMD_path);
end

% Find match of all datasets (default: day n and day n+1, can be changed to track between non-consecutive datasets)
for id = 1:numData-1
    input.data_path1 = ['D',num2str(id)]; % frist dataset, NEED CHANGE 
    input.data_path2 = ['D',num2str(id+1)]; % second dataset, NEED CHANGE 
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



%% ----------Plot matched units----------
% unit locations of two datasets, INPUT = result path, EMD file
plot_unit(chan_pos,input,'result12','EMD_post1.mat');

% waveform, INPUT = cluster number in dataset 1 and 2
plot_waveform(input,chan_pos,6,6);

% Z distance distribution
plot_z(input.result_path);

% waveform vs physical distance. INPUT = result path 
plot_dist(input,'result12');



%% Find chains
for id = 1:numData-1 % Load data
    all_input(id) = load(fullfile(input.input_path,['result',num2str(id),num2str(id+1)], "Input.mat")); 
    all_output(id) = load(fullfile(input.input_path,['result',num2str(id),num2str(id+1)],'Output.mat'));
end
[chain_all,z_loc,len] = chain_summary(all_input,all_output,numData);
save(fullfile(input.input_path,'chain_summary.mat'),"all_input","all_output",'chain_all','z_loc','len')
fprintf('Chains found. \n')



%----------Plot chains of interest (waveform, firing rate, original location, drift-corrected location, L2)----------
full_chain = chain_all(len == numData,:); %find chains with length across all datasets
[L2_weight,fr_all,fr_change,x_loc_all,z_loc_all] = chain_stats(all_input,all_output,full_chain,numData);

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


