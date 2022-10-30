from pathlib import Path
import tqdm

import numpy as np

from ibllib.ephys import neuropixel
from ibllib.io import spikeglx
import localization_pipeline.localizer
from localization_pipeline.localizer import LOCALIZER

REPO_PATH = Path(localization_pipeline.localizer.__file__).parents[1]
pids = ['8ca1a850-26ef-42be-8b28-c2e2d12f06d6', '8413c5c6-b42b-4ec6-b751-881a54413628', 'ce24bbe9-ae70-4659-9e9c-564d1a865de8', 'ce397420-3cd2-4a55-8fd1-5e28321981f4']
pid = pids[2]

standardized_file = next(Path(f"/datadisk/Data/spike_sorting/benchmark/raw/{pid}/").glob("*.ap.normalized.bin"))

h = neuropixel.trace_header(version=1)
residual_file = standardized_file
sr = spikeglx.Reader(str(standardized_file).replace('ap.normalized.bin', 'ap.cbin'))

params = dict(
    bin_file=standardized_file,
    # residual_file=standardized_file,
    dtype='float32',
    spike_index_path=standardized_file.parent.joinpath('detection', 'spike_index.npy'),
    templates_path=None,
    geom_path=np.c_[h['x'], h['y']],
    denoiser_weights=REPO_PATH.joinpath('pretrained_denoiser/denoise.pt'),
    offset_detector_denoiser=42,
    n_filters=[16, 8, 4],
    filter_sizes=[5, 11, 21],
    sampling_rate=sr.fs,
    multi_processing=1,
    n_processors=5,
    spike_size=121,
    n_channels_loc=10
)

localizer_obj = LOCALIZER(**params)

localizer_obj.load_denoiser()
output_directory = standardized_file.parent.joinpath('localisation')
if not output_directory.exists():
    output_directory.mkdir(exist_ok=True)
    for i in tqdm.tqdm(range(int(np.ceil(sr.rl)))):
        _ = localizer_obj.get_estimate(i, threshold=6, output_directory=output_directory)
