def sequences_to_fasta(filepath, sequences):
    with open(filepath, "w") as outfile:
        for i, sequence in enumerate(sequences):
            outfile.write(f">{i}\n")
            outfile.write(f"{sequence}\n")


def read_fasta(filepath):
    with open(filepath) as infile:
        next(infile)
        sequence = infile.read()
    return sequence.replace("\n", "")
