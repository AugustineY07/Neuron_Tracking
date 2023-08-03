# Neuron_Tracking
This repositor contains the code used for the manuscript: "Multi-day Neuron Tracking in High Density Electrophysiology Recordings using EMD".
All features are implemented in MATLAB R20221b.

Directories:
- Pipeline: code for tracking neurons
- Example: example code with 5 datasets from animal AL032
- Figure: code for reproducing figures in the manuscript

1. Pipeline
NT_main.m
Goal: 
Input: Kilosort cluster label, channel map, mean waveforms(preprocessed with ecephys_spike_sorting repo)  
Output: Unit match assignment  
Documentation:  
output.all_results column code = [reference, d2 clu label, d1 clu label, EMD distance, location distance, waveform distance, vertical distance]  










