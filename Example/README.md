An example with AL032 shank 1 data used in the manuscript.  

To run the example, please follow the steps below:  
  1. Clone the repo, or download a zip folder
  2. Add the parent folder (.../Neuron_Tracking) to the MATLAB folder with subfolders
  3. Edit paths and required parameters in 'Example_run.m' and run.  

Input required to run the pipleline with your own data:  
- channel_map.npy: standard KS2.5 output
- channel_positions.npy: standard KS2.5 output  
- cluster_KSLabel.tsv: standard KS2.5 output  
- ksproc_mean_waveforms.npy: mean waveforms in (nUnit x nChan x nTimepoints), saved as npy. In our examples we used C_Waves to calculate the mean waveforms. 
- metrics.csv: unit metrics; can be generated with ecephys_pipeline or other tool; only the first two columns (cluster ID and firing rate) are required. Note that cluster ID in the metrics table is 0-based; elsewhere cluster ID is zero based.

