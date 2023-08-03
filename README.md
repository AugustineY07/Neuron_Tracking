# Neuron_Tracking
This repositor contains the code used for the manuscript: "Multi-day Neuron Tracking in High Density Electrophysiology Recordings using EMD".
All features are implemented in MATLAB R20221b.

Directories:
- Pipeline: code for tracking neurons
- Example: example code with 5 datasets from animal AL032
- Figure: code for reproducing figures in the manuscript

1. Pipeline
Function 1 - NT_main.m  
Goal: Match units between 2 datasets.    
Input: Kilosort cluster label(only include clusters in the region of interest), channel map, mean waveforms(preprocessed with ecephys_spike_sorting repo)  
Output: Unit match assignment

User defined parameters:  
  
Documentation:  
output.all_results (column code) = [reference(0 for datasets without validation info), d2 clu label, d1 clu label, EMD distance, location distance, waveform distance, vertical distance]
output.results_wth = unit assignment with threshold applied  
  

    
2. Example
Function - Example_run.m  
Goal:
Input: Animal AL032 shank 1 dataset 1 to 5  
Output: 
  
  
  
3. Figure - code for reproducing figures 

