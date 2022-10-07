# Neuron_Tracking
This repositor keeps track of working progress of the project: "Neuron Tracking with Chronic NP2.0 Recordings from Mouse Visual Cortex".

Goal: 
Match neurons in chronic NP2.0 recordings after spike sorting using cluster location, spiking amplitude, template similarity and cluster spread difference. 

Data: 
1. Simulation data 
2. Neuropixels 2.0 Elctrophysiology Recording 

Method: 
3D Localization + EMD

Details:
The folder "mixed" contains simulated data with simutaneous adding and dropping points
computeSTAT.m calculates averaged template waveforms  
EMD_main.m generates and/or find matching pairs with simulation data using EMD algorithm
    Options - EMD package: matlab, CVX, TFOCS
              Dimension: 2D or 3D points
              Error of spike cluster estimation: 'Pos' = error in squared space, 'xyz' = error in a rectangle space
              Mode: how points are gegerated, 'standard' = uniform, 'lumpiness' = density differ by y-location, 'gainloss' = adding and/or dropping points 
