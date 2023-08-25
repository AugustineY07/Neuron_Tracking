% This is the main function of the algorithm
% Input: kilosort cluster label, channel map, mean waveforms, 
% Output: Unit match assignment
% For more comparisons, users need to write their own loops; see
% Example_run.m
function NT_main(input,chan_pos,mwf1,mwf2)

% Estimate location 
wf_metrics1 = wave_metrics(mwf1, chan_pos, input); %col 9,10,11 = x,z,y
wf_metrics2 = wave_metrics(mwf2, chan_pos, input);

% Estimate drift
output.threshold = input.threshold;
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
end