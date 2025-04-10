import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
from Bio import SeqIO
from tabulate import tabulate

# Read DNA sequence from FASTA file
file_path = "../data/sequence.fasta"
record = SeqIO.read(file_path, "fasta")
dna_sequence = str(record.seq)

# Nucleotide frequency
nucleotide_counts = Counter(dna_sequence)
df_frequencies = pd.DataFrame(nucleotide_counts.items(), columns=["Nucleotide", "Frequency"])

print(df_frequencies)

plt.figure(figsize=(6, 4))
sns.barplot(x="Nucleotide", y="Frequency", data=df_frequencies, palette="viridis")
plt.title("Nucleotide Frequency in DNA Sequence")
plt.xlabel("Nucleotide")
plt.ylabel("Count")
plt.show()

# DNA to mRNA (transcription)
def transcribe_dna_to_rna(dna_sequence):
    return dna_sequence.replace("T", "U")

rna_sequence = transcribe_dna_to_rna(dna_sequence)
print("mRNA sequence (first 100 bases):", rna_sequence[:100], "...")

# Codon table
codon_table = {
    'UUU': 'Phenylalanine', 'UUC': 'Phenylalanine', 'UUA': 'Leucine', 'UUG': 'Leucine',
    'UCU': 'Serine', 'UCC': 'Serine', 'UCA': 'Serine', 'UCG': 'Serine',
    'UAU': 'Tyrosine', 'UAC': 'Tyrosine', 'UAA': 'STOP', 'UAG': 'STOP', 'UGA': 'STOP',
    'UGU': 'Cysteine', 'UGC': 'Cysteine', 'UGG': 'Tryptophan',
    'CUU': 'Leucine', 'CUC': 'Leucine', 'CUA': 'Leucine', 'CUG': 'Leucine',
    'CCU': 'Proline', 'CCC': 'Proline', 'CCA': 'Proline', 'CCG': 'Proline',
    'CAU': 'Histidine', 'CAC': 'Histidine', 'CAA': 'Glutamine', 'CAG': 'Glutamine',
    'CGU': 'Arginine', 'CGC': 'Arginine', 'CGA': 'Arginine', 'CGG': 'Arginine',
    'AUU': 'Isoleucine', 'AUC': 'Isoleucine', 'AUA': 'Isoleucine', 'AUG': 'Methionine',  # START codon
    'ACU': 'Threonine', 'ACC': 'Threonine', 'ACA': 'Threonine', 'ACG': 'Threonine',
    'AAU': 'Asparagine', 'AAC': 'Asparagine', 'AAA': 'Lysine', 'AAG': 'Lysine',
    'AGU': 'Serine', 'AGC': 'Serine', 'AGA': 'Arginine', 'AGG': 'Arginine',
    'GUU': 'Valine', 'GUC': 'Valine', 'GUA': 'Valine', 'GUG': 'Valine',
    'GCU': 'Alanine', 'GCC': 'Alanine', 'GCA': 'Alanine', 'GCG': 'Alanine',
    'GAU': 'Aspartic acid', 'GAC': 'Aspartic acid', 'GAA': 'Glutamic acid', 'GAG': 'Glutamic acid',
    'GGU': 'Glycine', 'GGC': 'Glycine', 'GGA': 'Glycine', 'GGG': 'Glycine'
}

# Translation: mRNA to protein
def translate_rna_to_protein(rna_sequence):
    protein = []
    start_index = rna_sequence.find("AUG")
    if start_index == -1:
        return "No START codon (AUG) found in mRNA sequence."
    
    for i in range(start_index, len(rna_sequence) - 2, 3):
        codon = rna_sequence[i:i+3]
        aminoacid = codon_table.get(codon, "?")
        if aminoacid == "STOP":
            break
        protein.append(aminoacid)
    return "-".join(protein)

protein_sequence = translate_rna_to_protein(rna_sequence)
amino_counts = Counter(protein_sequence.split("-"))
df_amino = pd.DataFrame(amino_counts.items(), columns=["Amino Acid", "Frequency"])
df_amino = df_amino.sort_values(by="Frequency", ascending=False)

print(tabulate(df_amino, headers='keys', tablefmt='pretty'))

plt.figure(figsize=(10, 5))
sns.barplot(x="Amino Acid", y="Frequency", data=df_amino, palette="viridis")
plt.title("Amino Acid Frequency in Translated Protein")
plt.xlabel("Amino Acid")
plt.ylabel("Count")
plt.xticks(rotation=45)
plt.show()
