import os
import numpy as np

kilo_sort_dir = r"D:\Kilosort\Mouse08_SC_20251007_810to2250\kilosort4"
spike_times = np.load(os.path.join(kilo_sort_dir, "spike_times.npy"))
spike_clusters = np.load(os.path.join(kilo_sort_dir, "spike_clusters.npy"))

