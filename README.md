# Neuron_Tracking [Still working, please check back in a week]
This repositor contains the code used for the manuscript: "Multi-day Neuron Tracking in High Density Electrophysiology Recordings using EMD".
All features are implemented in MATLAB R2021b.

Required packages:  
1. npy-matlab(https://github.com/kwikteam/npy-matlab): need for accessing .npy files  
2. JRClust(https://github.com/jenniferColonell/JRCLUST): need for sorted dataset post-processing  

Directories:
- Pipeline: code for tracking neurons
- Example: example code with 5 datasets from animal AL032
- Figure: code for reproducing figures in the manuscript

1. Pipeline  

Function 1 - NT_main.m  
Goal: Match units between 2 datasets.    
Input: Kilosort cluster label(only include clusters in the region of interest), channel map, mean waveforms(preprocessed with ecephys_spike_sorting repo)  
Output: Unit match assignment  

Function 2 - chain_summary.mat  
Goal: Find chains using output from all datasets.    
Input: Unit assignment from multiple datasets   
Output: Summary of chains within distance threshold  

Function 3 - chain_stats.mat  
Goal: Summarize L2, FR, and locations of all fully tracked chains.   
Input: Summary of chains 
Output: Plots(waveforms, firing rate, xz locations) of selected chains  
  
Documentation:  
output.all_results (column code) = [reference(0 for datasets without validation info), d2 clu label, d1 clu label, EMD distance, location distance, waveform distance, vertical distance]
output.results_wth = unit assignment with threshold applied  
  

    
2. Example
Function - Example_run.m  
Goal: Demonstrate an example with 5 datasets
Input: Animal AL032 shank 1 dataset 1 to 5  
Output:
result_chain.mat - A matlab file contains the summary of chains within distance threshold  
chain_stats.mat - A matlab file contains the summary of L2 weight, firing rate, and xz locations of fully tracked chains  
Plot waveforms, firing rate, and xz locations of a chain of choice  
  
  
  
4. Figure - code for reproducing figures 

