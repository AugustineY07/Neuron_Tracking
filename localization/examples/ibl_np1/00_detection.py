from pathlib import Path
import shutil

import numpy as np

from ibllib.dsp.voltage import decompress_destripe_cbin
from ibllib.ephys import neuropixel
from ibllib.io import spikeglx
from ibllib.dsp.voltage import destripe
from ibllib.dsp.utils import rms

import detect.detector
from detect.run import run

APPLY_NN = False
BATCH_SIZE_SECS = 1  # for detection
REPO_PATH = Path(detect.detector.__file__).parents[1]
pids = ['8ca1a850-26ef-42be-8b28-c2e2d12f06d6', '8413c5c6-b42b-4ec6-b751-881a54413628', 'ce24bbe9-ae70-4659-9e9c-564d1a865de8', 'ce397420-3cd2-4a55-8fd1-5e28321981f4']
pid = pids[3]
cbin_file = next(Path(f"/datadisk/Data/spike_sorting/benchmark/raw/{pid}/").glob("*.ap.cbin"))

h = neuropixel.trace_header(version=1)
sr = spikeglx.Reader(cbin_file)

standardized_file = cbin_file.parent.joinpath(f"{cbin_file.stem}.normalized.bin")  # accounts fro both bin s and
if not standardized_file.exists():
    batch_size_secs = 1
    batch_intervals_secs = 50
    # scans the file at constant interval, with a demi batch starting offset
    nbatches = int(np.floor((sr.rl - batch_size_secs) / batch_intervals_secs - 0.5))
    wrots = np.zeros((nbatches, sr.nc - sr.nsync, sr.nc - sr.nsync))
    for ibatch in np.arange(nbatches):
        ifirst = int((ibatch + 0.5) * batch_intervals_secs * sr.fs + batch_intervals_secs)
        ilast = ifirst + int(batch_size_secs * sr.fs)
        sample = destripe(sr[ifirst:ilast, :-sr.nsync].T, fs=sr.fs, neuropixel_version=1)
        np.fill_diagonal(wrots[ibatch, :, :], 1 / rms(sample))
    wrot = np.median(wrots, axis=0)
    decompress_destripe_cbin(sr, h=h, wrot=wrot, output_file=standardized_file, dtype=np.float32, nc_out=sr.nc - sr.nsync)
    # shutil.copyfile(cbin_file.with_suffix('.meta'), standardized_file.with_suffix('.meta'))

params = dict(
    apply_nn=APPLY_NN,  # If set to false, run voltage threshold instead of NN detector,
    detect_threshold= .55 if APPLY_NN else 6,  # 0.5 if apply NN, 4/5/6 otherwise,
    filter_sizes_denoise=[5, 11, 21],
    geom_array=np.c_[h['x'], h['y']],
    len_recording=sr.rl,
    n_batches=sr.rl / 2,
    n_filters_denoise=[16, 8, 4],
    n_filters_detect=[16, 8, 8],
    n_processors=1,
    n_sec_chunk=BATCH_SIZE_SECS,
    n_sec_chunk_gpu_detect=.1,
    output_directory=cbin_file.parent.joinpath("detection"),
    path_nn_denoiser=REPO_PATH.joinpath('pretrained_denoiser/denoise.pt'),
    path_nn_detector=REPO_PATH.joinpath('pretrained_detector/detect_np1.pt'),
    run_chunk_sec='full',
    sampling_rate=sr.fs,
    spatial_radius=70,
    spike_size_nn=121,
    standardized_dtype='float32',
    standardized_path=standardized_file,
)

run(**params)
