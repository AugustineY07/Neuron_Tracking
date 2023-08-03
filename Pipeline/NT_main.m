% This is the main function of the algorithm
% Input: kilosort cluster label, channel map, mean waveforms, 
% Output: Unit match assignment
% For more comparisons, users need to write their own loops

addpath(genpath('C:\Users\labadmin\Desktop\Neuron Tracking Pipeline\User version'))
addpath(genpath('D:\Data\Pipeline\npy\npy-matlab'))

input.validation = 0; %whether there is validation dataset
input.validation_path = '';

% Input data
input.input_path = 'C:\Users\labadmin\Desktop\Neuron Tracking Pipeline\User version\'; %main directory 
input.EMD_path = fullfile(input.input_path,'EMD_input\'); %EMD input directory 
input.result_path = fullfile(input.input_path,'result1\'); %result directory 
input.fs = 30000; %acquisition rate
input.ts = 82; %wf time samples
input.l2_weights = 1500;
input.threshold = 10;
input.dim_mask = logical([1,1,1,0,0,0,0,0,0,1]);
input.dim_mask_physical = logical([1,1,1,0,0,0,0,0,0,0]);
input.dim_mask_wf = logical([0,0,0,0,0,0,0,0,0,1]);
input.data_path1 = 'D1';
input.data_path2 = 'D2';
input.chan_pos_name = 'channel_positions.npy';
input.wf_name = 'ksproc_mean_waveforms.npy';
input.KSLabel_name = 'cluster_KSLabel.tsv';

input.input_name = 'input1.mat';
input.input_name_post = 'input_post1.mat';
input.filename_pre = 'EMD_pre1.mat';
input.filename_post = 'EMD_post1.mat';

%ksLabel1 = readKS2label(fullfile(input_path,'D1\cluster_KSLabel.tsv'));
chan_pos = readNPY(fullfile(input.input_path,'D1\channel_positions.npy'));
mwf1 = readNPY(fullfile(input.input_path, 'D1\ksproc_mean_waveforms.npy'));
%ksLabel2 = readKS2label(fullfile(input_path,'D2\cluster_KSLabel.tsv'));
mwf2 = readNPY(fullfile(input.input_path, 'D2\ksproc_mean_waveforms.npy'));

% Estimate location 
wf_metrics1 = wave_metrics(mwf1, chan_pos, input); %col 9,10,11 = x,z,y
wf_metrics2 = wave_metrics(mwf2, chan_pos, input);

% Estimate drift
output.z_mode = 0;
output = create_EMD_input(input, output, wf_metrics1, wf_metrics2, mwf1, mwf2, 'pre'); 
output = EMD_unit_match(input,output,'pre');
[output.diffZ,edges] = z_estimate(input);
output.z_mode = kernelModeEstimate(output.diffZ);

% EMD unit matching
output = create_EMD_input(input, output, wf_metrics1, wf_metrics2, mwf1, mwf2, 'post'); 
output = EMD_unit_match(input,output,'post');

% thresholding
output.results_wth = output.all_results_post(output.all_results_post(:,7) <= input.threshold,:); 

% Save
if ~exist(input.result_path, 'dir')
    mkdir(input.result_path);
end
save(fullfile(input.result_path,'Input.mat'),"input");
save(fullfile(input.result_path,'Output.mat'),"output");
