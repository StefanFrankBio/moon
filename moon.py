#%%
import re
import collections
import statistics
import pandas as pd
import matplotlib.pyplot as plt


def read_fasta(filepath):
    with open(filepath) as infile:
        next(infile)
        sequence = infile.read()
    return sequence.replace("\n", "")


def read_multifasta(filepath):
    with open(filepath) as infile:
        peptides = infile.readlines()
    return [p.rstrip() for p in peptides if p[0] != ">"]


def findall_indices(regex, sequence):
    return [match.start() for match in re.finditer(regex, sequence)]


def split_by_modulo(int_list, modulo):
    split_list = [[] for _ in range(modulo)]
    for i in int_list:
        split_list[i % modulo].append(i)
    return split_list


def translate_dna(sequence):
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
    codons = re.findall("...", sequence)
    return "".join([codon_dict[c] for c in codons])


def findall(regex, sequence, overlap=False):
    if overlap:
        return re.findall(f"(?=({regex}))", sequence)
    else:
        return re.findall(f"{regex}", sequence)


def sequences_to_fasta(filepath, sequences):
    with open(filepath, "w") as outfile:
        for i, sequence in enumerate(sequences):
            outfile.write(f">{i}\n")
            outfile.write(f"{sequence}\n")


def relative_frequencies(outcome_list):
    number_of_outcomes = len(outcome_list)
    outcome_counter = collections.Counter(outcome_list)
    return {k: v / number_of_outcomes for k, v in outcome_counter.items()}


def mean_frequency(frequency_dict, key_list):
    frequencies = [frequency_dict[key] for key in key_list]
    return statistics.mean(frequencies)


def reverse_complement(sequence):
    trans_dict = str.maketrans("ATGC", "TACG")
    return sequence.translate(trans_dict)[::-1]


def sequence_neighborhood(sequence, position, offset):
    if offset > position:
        offset = position
    return sequence[position - offset : position + offset]


def flatten(list_of_sublists):
    return [item for sublist in list_of_sublists for item in sublist]


def ratio_dict(dict1, dict2):
    return {k: dict1[k] / dict2[k] for k in dict1.keys()}


def sort_dict_by_value(dictionary, reverse=False):
    return dict(sorted(dictionary.items(), key=lambda item: item[1], reverse=reverse))


def build_orf_table(sequence):
    stop_codon_positions = findall_indices("TAA|TAG|TGA", sequence)
    stop_codon_positions_by_frame = split_by_modulo(stop_codon_positions, 3)
    orf_table = pd.DataFrame()
    orf_table["upstream_end"] = [
        i + 3 for frame in stop_codon_positions_by_frame for i in frame[:-1]
    ]
    orf_table["downstream_end"] = [
        i for frame in stop_codon_positions_by_frame for i in frame[1:]
    ]
    orf_table["dna_sequence"] = [
        sequence[i:j] for i, j in zip(orf_table.upstream_end, orf_table.downstream_end)
    ]
    orf_table["amino_sequence"] = [translate_dna(dna) for dna in orf_table.dna_sequence]
    orf_table["length"] = list(map(len, orf_table.dna_sequence))
    orf_mean_length = int(orf_table.length.mean())  # mean_length = 39
    orf_table = orf_table[orf_table.length > orf_mean_length]
    codons = findall("...", sequence, overlap=True)
    relative_codon_frequencies = relative_frequencies(codons)
    codons_by_orf = [findall("...", seq) for seq in orf_table.dna_sequence]
    orf_table["mean_codon_frequency"] = [
        mean_frequency(relative_codon_frequencies, codons) for codons in codons_by_orf
    ]
    orf_table["upstream_end_neighborhood_sequence"] = [
        sequence_neighborhood(sequence, i, orf_mean_length)
        for i in orf_table.upstream_end
    ]
    return orf_table


def kmer_analysis(sequence, orf_table, kmer_length):
    kmers = findall("." * kmer_length, sequence, overlap=True)
    kmer_counter = relative_frequencies(kmers)
    neighborhood_kmers = [
        findall("." * kmer_length, n, overlap=True)
        for n in orf_table[orf_table.length >= 1000].upstream_end_neighborhood_sequence
    ]
    neighborhood_kmers = flatten(neighborhood_kmers)
    neighborhood_kmer_counter = relative_frequencies(neighborhood_kmers)
    kmer_ratio_dict = ratio_dict(neighborhood_kmer_counter, kmer_counter)
    kmer_ratio_dict = sort_dict_by_value(kmer_ratio_dict, reverse=True)
    neighborhood_kmers_total = [
        findall("." * kmer_length, n, overlap=True)
        for n in orf_table.upstream_end_neighborhood_sequence
    ]
    orf_table["max_kmer_ratio"] = [
        max([kmer_ratio_dict.get(key, 0) for key in kmer_list])
        for kmer_list in neighborhood_kmers_total
    ]
    return orf_table


if __name__ == "__main__":
    sequence = read_fasta("data/NC_009641.fasta")
    orf_table = build_orf_table(sequence)
    sequence_complement = reverse_complement(sequence)
    orf_table_complement = build_orf_table(sequence_complement)
    orf_table = kmer_analysis(sequence, orf_table, 5)
    orf_table_complement = kmer_analysis(sequence_complement, orf_table_complement, 5)
    orf_table = pd.concat([orf_table, orf_table_complement])
    orf_table.to_csv("data/orf_table.tsv", sep="\t", index=False)

    peptides = read_multifasta("data/filtered_peptides.fasta")
    checklist = [0 for _ in range(len(orf_table.amino_sequence))]
    for i, seq in enumerate(orf_table.amino_sequence):
        if any([pep in seq for pep in peptides]):
            checklist[i] = 1
    orf_table["peptide_found"] = checklist
    orf_table.to_csv("data/orf_table.tsv", sep="\t", index=False)
