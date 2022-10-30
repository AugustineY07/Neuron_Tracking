#!/usr/bin/env python
# coding: utf-8

import logging
import os
try:
    from pathlib2 import Path
except ImportError:
    from pathlib import Path

import numpy as np
import torch
import torch.multiprocessing as mp
import parmap
#import spikes_localization_registration 

from detect.detector import Detect
from localization_pipeline.denoiser import Denoise

from detect.deduplication import deduplicate_gpu, deduplicate

from scipy.signal import argrelmin

from detect.run import run

# In[22]: 
    
    
# parameters and input path
geom_path = 'Z:/AL032_out/results/catgt_AL032_2019-12-03_stripe192_natIm_sh_g0/AL032_2019-12-03_stripe192_natIm_sh_g0_imec0/imec0_ks25/channel_positions.npy'
path_nn_detector = 'D:/Data/Papers/3D Localization/spikes_localization_registration/pretrained_detector/detect_np2.pt'
path_nn_denoiser = 'D:/Data/Papers/3D Localization/spikes_localization_registration/pretrained_denoiser/denoise.pt'
standardized_path = 'Z:/AL032_out/results/catgt_AL032_2019-12-03_stripe192_natIm_sh_g0/AL032_2019-12-03_stripe192_natIm_sh_g0_imec0/Localization/standardized.bin'
standardized_dtype = 'float32'
sampling_rate = 30000
len_recording = 1917 #1964/1917/2420/2273
output_directory = 'Z:/AL032_out/results/catgt_AL032_2019-12-03_stripe192_natIm_sh_g0/AL032_2019-12-03_stripe192_natIm_sh_g0_imec0/Localization'

geom_array = np.load(geom_path)
apply_nn = False ### If set to false, run voltage threshold instead of NN detector 
spatial_radius = 70 #what is this?
n_sec_chunk = 1
n_processors = 1
n_sec_chunk_gpu_detect = .1
detect_threshold = 4 ## 0.5 if apply NN, 4/5/6 otherwise 
n_filters_detect = [16, 8, 8] 
spike_size_nn = 121 ### In sample steps, not used 
n_filters_denoise = [16, 8, 4]
filter_sizes_denoise = [5, 11, 21]
n_batches = len_recording//n_sec_chunk




#print(localizer.read_waveforms)

# In[ ]:


sampling_rate = 30000

run(standardized_path, standardized_dtype, output_directory, geom_array, spatial_radius, apply_nn, n_sec_chunk, n_batches, n_processors, n_sec_chunk_gpu_detect, sampling_rate, len_recording,
    detect_threshold, path_nn_detector, n_filters_detect, spike_size_nn, path_nn_denoiser, n_filters_denoise, filter_sizes_denoise, run_chunk_sec='full')


# In[ ]:


#### LOCALIZATION

import os
import numpy as np
from tqdm import tqdm
# from residual import RESIDUAL
from localization_pipeline.localizer import LOCALIZER
from localization_pipeline.merge_results import get_merged_arrays


# In[ ]:


### Change paths to your data here
bin_file = standardized_path
residual_file = bin_file
dtype_input = 'float32'

fname_spike_train = 'Z:/AL032_out/results/catgt_AL032_2019-12-03_stripe192_natIm_sh_g0/AL032_2019-12-03_stripe192_natIm_sh_g0_imec0//Localization/spike_channel.npy'
# Sort spike train if not 
spt_array = np.load(fname_spike_train)
spt_array = spt_array.astype('int')
spt_array = spt_array[spt_array[:, 0].argsort()]
np.save(fname_spike_train, spt_array)


#geom_path = "D:/Pipeline/Localization/day3 all batch/geom.npy"
n_channels = np.load(geom_path).shape[0]

denoiser_weights = path_nn_denoiser
denoiser_min = 42 ## Goes with the weights

#n_batches = 1000
#len_recording = 120 #1964
#sampling_rate = 30000
n_channels_loc = 10

fname_templates = None


# In[ ]:

                         #standardized  float32  spike-and-max-channel  none         channel-map-pretrained-denoiser peak-alignment
localizer_obj = LOCALIZER(bin_file, dtype_input, fname_spike_train, fname_templates, geom_path, denoiser_weights, denoiser_min)
#localizer_obj = LOCALIZER(bin_file, dtype_input, fname_spike_train, fname_templates, geom_path, denoiser_weights, denoiser_min, n_filters_denoise, filter_sizes_denoise, sampling_rate, spike_size_nn, n_channels_loc)
# localizer_obj.get_offsets()
# localizer_obj.compute_aligned_templates()
localizer_obj.load_denoiser()
output_directory = 'Z:/AL032_out/results/catgt_AL032_2019-12-03_stripe192_natIm_sh_g0/AL032_2019-12-03_stripe192_natIm_sh_g0_imec0//Localization/position_results_files'
if not os.path.exists(output_directory):
    os.makedirs(output_directory)
merge_directory = 'Z:/AL032_out/results/catgt_AL032_2019-12-03_stripe192_natIm_sh_g0/AL032_2019-12-03_stripe192_natIm_sh_g0_imec0//Localization/position_results_files_merged'      
if not os.path.exists(merge_directory):
    os.makedirs(merge_directory)
for i in tqdm(range(n_batches)):
    localizer_obj.get_estimate(i, threshold = detect_threshold, output_directory = 'Z:/AL032_out/results/catgt_AL032_2019-12-03_stripe192_natIm_sh_g0/AL032_2019-12-03_stripe192_natIm_sh_g0_imec0//Localization/position_results_files')


# In[10]:


from localization_pipeline.merge_results import get_merged_arrays


# In[12]:

get_merged_arrays('Z:/AL032_out/results/catgt_AL032_2019-12-03_stripe192_natIm_sh_g0/AL032_2019-12-03_stripe192_natIm_sh_g0_imec0//Localization/position_results_files', 'Z:/AL032_out/results/catgt_AL032_2019-12-03_stripe192_natIm_sh_g0/AL032_2019-12-03_stripe192_natIm_sh_g0_imec0//Localization/position_results_files_merged', n_batches)


# In[13]:

#from IPython import get_ipython

#get_ipython().run_line_magic('matplotlib', 'inline')

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec as gridspec
import matplotlib.cm as cm


# In[28]:

os.environ['KMP_DUPLICATE_LIB_OK']='True'
ptp_array = np.load('Z:/AL032_out/results/catgt_AL032_2019-12-03_stripe192_natIm_sh_g0/AL032_2019-12-03_stripe192_natIm_sh_g0_imec0//Localization/position_results_files_merged/results_max_ptp_merged.npy')
idx_good = np.where(ptp_array != 0)[0]
x_results = np.load('Z:/AL032_out/results/catgt_AL032_2019-12-03_stripe192_natIm_sh_g0/AL032_2019-12-03_stripe192_natIm_sh_g0_imec0//Localization/position_results_files_merged/results_x_merged.npy')[idx_good]
y_results = np.load('Z:/AL032_out/results/catgt_AL032_2019-12-03_stripe192_natIm_sh_g0/AL032_2019-12-03_stripe192_natIm_sh_g0_imec0//Localization/position_results_files_merged/results_y_merged.npy')[idx_good]
z_results = np.load('Z:/AL032_out/results/catgt_AL032_2019-12-03_stripe192_natIm_sh_g0/AL032_2019-12-03_stripe192_natIm_sh_g0_imec0//Localization/position_results_files_merged/results_z_merged.npy')[idx_good]
alpha_results = np.load('Z:/AL032_out/results/catgt_AL032_2019-12-03_stripe192_natIm_sh_g0/AL032_2019-12-03_stripe192_natIm_sh_g0_imec0//Localization/position_results_files_merged/results_alpha_merged.npy')[idx_good]
times_array = np.load('Z:/AL032_out/results/catgt_AL032_2019-12-03_stripe192_natIm_sh_g0/AL032_2019-12-03_stripe192_natIm_sh_g0_imec0//Localization/position_results_files_merged/times_read.npy')[idx_good]
channelMap = np.load('Z:/AL032_out/results/catgt_AL032_2019-12-03_stripe192-natIm_sh_g0/AL032_2019-12-03_stripe192-natIm_sh_g0_imec0/imec0_ks25/channel_positions.npy')


# In[21]:


plt.hist(ptp_array[idx_good], bins = 50)
plt.show()


# In[24]:


ptp_rescaled = ptp_array[idx_good] - ptp_array[idx_good].min()
ptp_rescaled = ptp_rescaled/ptp_rescaled.max()
ptp_rescaled[ptp_rescaled >= 0.4] = 0.4
ptp_rescaled = ptp_rescaled/ptp_rescaled.max() #make color seperation 

# In[25]:


plt.hist(ptp_rescaled, bins = 50)
plt.show()


# In[30]:
#plot enpty panal
fig = plt.figure(figsize = (26, 20))
spec = gridspec(ncols=8, nrows=1, figure=fig)

vir = cm.get_cmap('viridis')
f_ax = fig.add_subplot(spec[0, 1])
f_ax.scatter(channelMap[:, 0], channelMap[:, 1], c = 'orange', label = "NP channels", marker = 's', s = 10) #plot channels
f_ax.set_xlim((-100, 100))
    
# In[31]:


fig = plt.figure(figsize = (26, 20))
spec = gridspec(ncols=8, nrows=1, figure=fig)

vir = cm.get_cmap('viridis')
f_ax = fig.add_subplot(spec[0, 1])

#f_ax.set_xlim((-40, 62))
#f_ax.scatter(x_results, z_results, s = 2, color = vir(ptp_rescaled), alpha = 0.05) #vir(ptp_scaled_high_units) 0.05
f_ax.scatter(x_results, z_results, s = 2, color = vir(ptp_rescaled), alpha = 0.05)
f_ax.scatter(channelMap[:, 0], channelMap[:, 1], c = 'orange', label = "NP channels", marker = 's', s = 10) #plot channels

#f_ax.set_yticks([])
#f_ax.set_xticklabels(np.arange(-40, 62, 20), fontsize = 18)
#f_ax.set_ylim(2500, 3800)

f_ax.set_xlabel("X", fontsize = 22)

f_ax.set_title("LS position \n by max-ptp", fontsize = 18)
plt.show()


# In[ ]:

fig = plt.figure(figsize = (26, 20))
spec = gridspec(ncols=8, nrows=1, figure=fig)

vir = cm.get_cmap('viridis')
f_ax = fig.add_subplot(spec[0, 1])

#f_ax.set_xlim((-40, 62))
f_ax.scatter(y_results, z_results, s = 2, color = vir(ptp_rescaled), alpha = 0.05) #vir(ptp_scaled_high_units) 0.05
f_ax.scatter(channelMap[:, 0], channelMap[:, 1], c = 'orange', label = "NP channels", marker = 's', s = 10) #plot channels

#f_ax.set_yticks([])
#f_ax.set_xticklabels(np.arange(-40, 62, 20), fontsize = 18)
#f_ax.set_ylim(2500, 3800)

f_ax.set_xlabel("Y", fontsize = 22)

f_ax.set_title("LS position \n by max-ptp", fontsize = 18)
plt.show()


# In[ ]:

fig = plt.figure(figsize = (26, 20))
spec = gridspec(ncols=8, nrows=1, figure=fig)

vir = cm.get_cmap('viridis')
f_ax = fig.add_subplot(spec[0, 1])

#f_ax.set_xlim((-40, 62))
f_ax.scatter(alpha_results, z_results, s = 2, color = vir(ptp_rescaled), alpha = 0.05) #vir(ptp_scaled_high_units) 0.05
f_ax.scatter(channelMap[:, 0], channelMap[:, 1], c = 'orange', label = "NP channels", marker = 's', s = 10) #plot channels

#f_ax.set_yticks([])
#f_ax.set_xticklabels(np.arange(-40, 62, 20), fontsize = 18)
#f_ax.set_ylim(2500, 3800)

f_ax.set_xlabel("Alpha", fontsize = 22)

f_ax.set_title("LS position \n by max-ptp", fontsize = 18)
plt.show()


