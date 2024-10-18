# Neuron_Tracking  
## This repository has been moved to: https://github.com/janelia-TDHarrisLab/Yuan-Neuron_Tracking, please refer to it for the latest updates

This repository contains the code used for the manuscript: "Multi-day Neuron Tracking in High Density Electrophysiology Recordings using EMD". (doi: https://doi.org/10.1101/2023.08.03.551724)

All features are implemented in MATLAB R2023b.

## Required packages:  
* npy-matlab(https://github.com/kwikteam/npy-matlab): need for reading .npy files
* C_Waves (https://billkarsh.github.io/SpikeGLX/#post-processing-tools): to generate mean waveforms
* ecephys(https://github.com/jenniferColonell/ecephys_spike_sorting): for KS2.5 post-processing and metrics 

## Directories:
- Pipeline: code for tracking neurons
- Example: example code with 5 datasets from animal AL032
- Figure: code for reproducing figures in the manuscript  



### 1. Pipeline  

* Function 1 - NT_main.m  
Goal: Match units between 2 datasets.    
Input: Kilosort cluster label(only include clusters in the region of interest), channel map, mean waveforms (created with C_Waves or other tool, from the KS2.5 drift corrected data)  
Output: Unit match assignment  

* Function 2 - chain_summary.mat  
Goal: Find chains using output from all datasets.    
Input: Unit assignment from multiple datasets   
Output: Summary of chains within distance threshold  

* Function 3 - chain_stats.mat  
Goal: Summarize L2, FR, and locations of all fully tracked chains.   
Input: Summary of chains 
Output: Plots(waveforms, firing rate, xz locations) of selected chains  

* Additional Functions:
  * ks2_working_toBinary: Converts the KS2.5 temp_wh.dat file, which has been drift corrected, to a standard, correctly scaled binary to be input to C_Waves
  * dist_mat.m: A function used to generate a 'distance matrix' that gives an idea of whether to exclude certain datasets with poor trackability. It is more meaningful to run on datasets with many comparisons. 
  * acc.m: A function used to compute recovery rate(without threshold) and accuracy(with threshold). Only use if there are validation information. Will use a matrix with correct-match cluster labels as input. The directory 'supplementary' contains one example validation table 'truth.mat'(from animal AL032 shank 1 dataset 1 and 2).   


* Documentation:  
output.all_results (column code) = [reference(0 for datasets without validation info), d2 clu label, d1 clu label, EMD distance, location distance, waveform distance, vertical distance]
output.results_wth = unit assignment with threshold applied  
  

    
### 2. Example
Function - Example_run.m  
Goal: Demonstrate an example with 5 datasets
Input: Animal AL032 shank 1 dataset 1 to 5  
Output:
* result_chain.mat - A matlab file contains the summary of chains within distance threshold  
* chain_stats.mat - A matlab file contains the summary of L2 weight, firing rate, and xz locations of fully tracked chains  
* Plot waveforms, firing rate, and xz locations of a chain of choice  
  
  
  
### 3. Figure - code for reproducing figures 




