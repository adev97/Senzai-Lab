import os
import numpy as np
import matplotlib.pyplot as plt


kilo_sort_dir = r"D:\Kilosort\Mouse08_SC_20251007_810to2250\kilosort4"
os.chdir(kilo_sort_dir)

chanMap = np.load('channel_map.npy')
pos = np.load('channel_positions.npy')
# print(type(chanMap))
print(chanMap.shape)
print(pos.shape)
# print(chanMap.dtype)

## Visualize channel layout
x = pos[:, 0]
y = pos[:, 1]

plt.figure(figsize=(4,8))
plt.scatter(x, y, c=chanMap, s=20)
plt.gca().invert_yaxis()   # depth increases downward
plt.axis('equal')
plt.colorbar(label='Channel index (0-based)')
plt.xlabel('X (µm)')
plt.ylabel('Y (µm)')
plt.title('Kilosort4 Channel Map')
plt.show()

order = np.argsort(y)
print("Top of probe:")
print(chanMap[order][:20])

print("\nBottom of probe:")
print(chanMap[order][-20:])


# x = chmap['xcoords']
# y = chmap['ycoords']
# chan = chmap['chanMap']
# connected = chmap.get('connected', np.ones_like(chan, dtype=bool))

# plt.figure()
# plt.scatter(x[connected], y[connected], c=chan[connected], s=40)
# plt.gca().invert_yaxis()   # depth increases downward
# plt.axis('equal')
# plt.colorbar(label='Channel index')
# plt.xlabel('X (µm)')
# plt.ylabel('Y (µm)')
# plt.title('Kilosort Channel Map')
# plt.show()





























# load kilosort spike_times and spike_clusters
# kilo_sort_dir = r"D:\Kilosort\Mouse08_SC_20251007_810to2250\kilosort4"
# spike_times = np.load(os.path.join(kilo_sort_dir, "spike_times.npy"))
# spike_clusters = np.load(os.path.join(kilo_sort_dir, "spike_clusters.npy"))

# load continuous raw data


