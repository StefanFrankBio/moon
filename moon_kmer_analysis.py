import moontools
import re
import collections
import pandas as pd
import matplotlib.pyplot as plt

sequence = moontools.read_fasta("NC_009641.fasta")
kmers = re.findall(r"(?=(.....))", sequence)
kmer_counter = collections.Counter(kmers)
for key in kmer_counter.keys():
    kmer_counter[key] /= len(kmers)

df = pd.read_csv("orf_table.tsv", sep="\t")
stop_missing_sequences = list(df[df.length >= 1000].upstream_end_neighborhood_sequence)
stop_missing_kmers = []
for seq in stop_missing_sequences:
    if type(seq) != float:
        stop_missing_kmers += re.findall(r"(?=(.....))", seq)
stop_missing_kmer_counter = collections.Counter(stop_missing_kmers)
for key in stop_missing_kmer_counter.keys():
    stop_missing_kmer_counter[key] /= len(stop_missing_kmers)
delta_kmer_frequencies = {
    key: stop_missing_kmer_counter[key] / kmer_counter[key]
    for key in stop_missing_kmer_counter.keys()
}
delta_kmer_frequencies = dict(
    sorted(delta_kmer_frequencies.items(), key=lambda item: item[1], reverse=True)
)
kmer_positions = list(range(80 - 4)) * len(stop_missing_sequences)

kmer_dict = {key: [0 for _ in range(80)] for key in stop_missing_kmer_counter.keys()}

for i, kmer in enumerate(stop_missing_kmers):
    pos = kmer_positions[i]
    for j in range(5):
        kmer_dict[kmer][pos + j] += 1

plot_x = range(-40, 40)
plot_y_list = []
for k, v in delta_kmer_frequencies.items():
    if v > 5:
        plot_y_list.append(kmer_dict[k])
cum_list = [sum(x) for x in zip(*plot_y_list)]
plt.plot(plot_x, cum_list)
plt.show()