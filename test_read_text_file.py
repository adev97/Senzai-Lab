import numpy as np

text = np.load(r'D:\OpenEphys_Data\Mouse08\Mouse08_20251007_4shanks_810to2250_RFMapping\Record Node 102\experiment1\recording1\events\MessageCenter\text.npy')

print(type(text))
print(text.dtype)
print(text.shape)

sample_numbers = np.load(r'D:\OpenEphys_Data\Mouse08\Mouse08_20251007_4shanks_810to2250_RFMapping\Record Node 102\experiment1\recording1\events\MessageCenter\sample_numbers.npy')
print(type(sample_numbers))
print(sample_numbers.dtype)
print(sample_numbers.shape)

sample_numbers_from_events = np.load(r'D:\OpenEphys_Data\Mouse08\Mouse08_20251007_4shanks_810to2250_RFMapping\Record Node 102\experiment1\recording1\events\OneBox-100.OneBox-ADC\TTL\sample_numbers.npy')
print(type(sample_numbers_from_events))
print(sample_numbers_from_events.dtype)
print(sample_numbers_from_events.shape)

full_words = np.load(r"D:\OpenEphys_Data\Mouse08\Mouse08_20251007_4shanks_810to2250_RFMapping\Record Node 102\experiment1\recording1\events\OneBox-100.OneBox-ADC\TTL\full_words.npy")
