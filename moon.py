import re
import statistics
import collections
import pandas as pd


with open("NC_009641.fasta") as infile:
    next(infile)
    sequence = infile.read()
sequence = sequence.replace("\n", "")
stop_positions = [match.start() for match in re.finditer("TAA|TAG|TGA", sequence)]
stop_positions_by_frame = [[], [], []]
for i in stop_positions:
    stop_positions_by_frame[i % 3].append(i)
upstream_ends = [i + 3 for frame in stop_positions_by_frame for i in frame[:-1]]
downstream_ends = [i for frame in stop_positions_by_frame for i in frame[1:]]
dna_sequences = [sequence[i:j] for i, j in zip(upstream_ends, downstream_ends)]
codon_dict = {
    "AAA": "K",
    "AAC": "N",
    "AAG": "K",
    "AAT": "N",
    "ACA": "T",
    "ACC": "T",
    "ACG": "T",
    "ACT": "T",
    "AGA": "R",
    "AGC": "S",
    "AGG": "R",
    "AGT": "S",
    "ATA": "I",
    "ATC": "I",
    "ATG": "M",
    "ATT": "I",
    "CAA": "Q",
    "CAC": "H",
    "CAG": "Q",
    "CAT": "H",
    "CCA": "P",
    "CCC": "P",
    "CCG": "P",
    "CCT": "P",
    "CGA": "R",
    "CGC": "R",
    "CGG": "R",
    "CGT": "R",
    "CTA": "L",
    "CTC": "L",
    "CTG": "L",
    "CTT": "L",
    "GAA": "E",
    "GAC": "D",
    "GAG": "E",
    "GAT": "D",
    "GCA": "A",
    "GCC": "A",
    "GCG": "A",
    "GCT": "A",
    "GGA": "G",
    "GGC": "G",
    "GGG": "G",
    "GGT": "G",
    "GTA": "V",
    "GTC": "V",
    "GTG": "V",
    "GTT": "V",
    "TAA": "_",
    "TAC": "Y",
    "TAG": "_",
    "TAT": "Y",
    "TCA": "S",
    "TCC": "S",
    "TCG": "S",
    "TCT": "S",
    "TGA": "_",
    "TGC": "C",
    "TGG": "W",
    "TGT": "C",
    "TTA": "L",
    "TTC": "F",
    "TTG": "L",
    "TTT": "F",
}
amino_sequences = []
for seq in dna_sequences:
    codons = re.findall("...", seq)
    amino_sequences.append("".join([codon_dict[c] for c in codons]))
lengths = [j - i for i, j in zip(upstream_ends, downstream_ends)]
codons = []
for i in range(3):
    codons += re.findall("...", sequence[i:])
codon_counter = collections.Counter(codons)
len_codons = len(codons)
for key in codon_counter.keys():
    codon_counter[key] /= len(codons)
codon_frequency_total = statistics.mean([codon_counter[codon] for codon in codons])
codon_frequencies = []
for seq in dna_sequences:
    codons = re.findall("...", seq)
    codon_frequencies.append([codon_counter[codon] for codon in codons])
mean_codon_frequencies = [
    statistics.mean(freqs) if len(freqs) > 0 else 0 for freqs in codon_frequencies
]
stdev_codon_frequencies = [
    statistics.stdev(freqs) if len(freqs) > 1 else 0 for freqs in codon_frequencies
]

delta_from_mean_codon_frequencies = [
    abs(codon_frequency_total - i) for i in mean_codon_frequencies
]

upstream_end_neighborhood_sequences = [sequence[i - 40 : i + 40] for i in upstream_ends]

df = pd.DataFrame(
    {
        "upstream_end": upstream_ends,
        "downstream_end": downstream_ends,
        "upstream_end_neighborhood_sequence": upstream_end_neighborhood_sequences,
        "length": lengths,
        "mean_codon_frequency": mean_codon_frequencies,
        "stdev_codon_frequency": stdev_codon_frequencies,
        "delta_codon_frequency": delta_from_mean_codon_frequencies,
        "dna_sequence": dna_sequences,
        "amino_sequence": amino_sequences,
    }
)

df.to_csv("orf_table.tsv", sep="\t", index=False)

df = pd.read_csv("orf_table.tsv", sep="\t")
