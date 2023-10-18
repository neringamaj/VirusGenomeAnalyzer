from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
import math
import numpy as np


codons = {
    "ATA":"I", "ATC":"I", "ATT":"I", "ATG":"M",
    "ACA":"T", "ACC":"T", "ACG":"T", "ACT":"T",
    "AAC":"N", "AAT":"N", "AAA":"K", "AAG":"K",
    "AGC":"S", "AGT":"S", "AGA":"R", "AGG":"R",
    "CTA":"L", "CTC":"L", "CTG":"L", "CTT":"L",
    "CCA":"P", "CCC":"P", "CCG":"P", "CCT":"P",
    "CAC":"H", "CAT":"H", "CAA":"Q", "CAG":"Q",
    "CGA":"R", "CGC":"R", "CGG":"R", "CGT":"R",
    "GTA":"V", "GTC":"V", "GTG":"V", "GTT":"V",
    "GCA":"A", "GCC":"A", "GCG":"A", "GCT":"A",
    "GAC":"D", "GAT":"D", "GAA":"E", "GAG":"E",
    "GGA":"G", "GGC":"G", "GGG":"G", "GGT":"G",
    "TCA":"S", "TCC":"S", "TCG":"S", "TCT":"S",
    "TTC":"F", "TTT":"F", "TTA":"L", "TTG":"L",
    "TAC":"Y", "TAT":"Y", "TAA":"_", "TAG":"_",
    "TGC":"C", "TGT":"C", "TGA":"_", "TGG":"W"
}

START_CODON = "ATG"
STOP_CODONS = ["TAA", "TAG", "TGA"]

def get_start_stop_pairs(seq):
    pairs = []

    def find_pairs(dna_sequence):
        start_index = dna_sequence.find(START_CODON)
        while start_index != -1:
            seq_from_start = dna_sequence[start_index:]
            stop_index = -1
            for stop_codon in STOP_CODONS:
                idx = seq_from_start.find(stop_codon)
                if idx != -1:
                    if stop_index == -1 or idx < stop_index:
                        stop_index = idx

            if stop_index != -1:
                between_sequence = dna_sequence[start_index:start_index+stop_index]
                if all(between_sequence.find(stop) == -1 for stop in STOP_CODONS):
                    pairs.append((start_index, start_index + stop_index))

            start_index = dna_sequence.find(START_CODON, start_index + 1)

    find_pairs(seq)
    find_pairs(str(Seq(seq).reverse_complement()))

    return pairs

def compute_codon_dicodon_frequencies(sequence):
    codon_freq = defaultdict(int)
    dicodon_freq = defaultdict(int)
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i+3]
        if codon in codons:  # Ensure the codon exists in our dictionary
            codon_freq[codon] += 1
        if i + 6 <= len(sequence):
            dicodon = sequence[i:i+6]
            dicodon_freq[dicodon] += 1
    return codon_freq, dicodon_freq

# Determine variance across sequences
def variance_across_sequences(sequences, n):
    all_freqs = defaultdict(list)
    
    # Get frequencies for each sequence
    for seq in sequences:
        codon_freq, dicodon_freq = compute_codon_dicodon_frequencies(seq)
        if n == 3:
            freq = codon_freq
        else:
            freq = dicodon_freq
        
        for codon, count in freq.items():
            all_freqs[codon].append(count)
    
    # Calculate variance for each codon/dicodon
    variances = {}
    for codon, counts in all_freqs.items():
        variances[codon] = np.var(counts)

    # Sort by variance
    sorted_variances = sorted(variances.items(), key=lambda x: x[1], reverse=True)
    return sorted_variances

def compute_euclidean_distance(freq1, freq2):
    keys = set(freq1.keys()).union(set(freq2.keys()))
    distance = sum([(freq1[k] - freq2[k])**2 for k in keys])
    return math.sqrt(distance)

def compute_distance_matrix(sequences, use_dicodon=False):
    sequences = [seq for seq in sequences if len(seq) >= 100]
    
    precomputed_freqs = [compute_codon_dicodon_frequencies(seq) for seq in sequences]
    
    matrix = []
    for i, seq in enumerate(sequences):
        row = []
        for j, seq in enumerate(sequences):
            freq1 = precomputed_freqs[i][1 if use_dicodon else 0]
            freq2 = precomputed_freqs[j][1 if use_dicodon else 0]
            row.append(compute_euclidean_distance(freq1, freq2))
        matrix.append(row)
    return matrix

def format_phylip(matrix, names):
    n = len(matrix)
    output = [str(n)]
    for i in range(n):
        output.append(names[i] + " " + " ".join(map(str, matrix[i])))
    return "\n".join(output)

# Main execution
file_names = ["mamalian1.fasta", "mamalian2.fasta", "mamalian3.fasta", "mamalian4.fasta", 
              "bacterial1.fasta", "bacterial2.fasta", "bacterial3.fasta", "bacterial4.fasta"]

all_sequences = []
names = []
for file_name in file_names:
    for record in SeqIO.parse(file_name, "fasta"):
        names.append(record.id)
        seq = str(record.seq)
        protein_seq = []
        for start, stop in get_start_stop_pairs(seq):
            protein_seq.extend([codons[seq[i:i+3]] for i in range(start, stop, 3)])
        all_sequences.append("".join(protein_seq))

codon_distance_matrix = compute_distance_matrix(all_sequences)
dicodon_distance_matrix = compute_distance_matrix(all_sequences, use_dicodon=True)


codon_variances = variance_across_sequences(all_sequences, 3)
dicodon_variances = variance_across_sequences(all_sequences, 6)

print("Codons with the most variation:", codon_variances[:5])
print("Dicodons with the most variation:", dicodon_variances[:5])

codon_phylip_output = format_phylip(codon_distance_matrix, names)
dicodon_phylip_output = format_phylip(dicodon_distance_matrix, names)

with open('codon_distance_matrix.txt', 'w') as f:
    f.write(codon_phylip_output)

with open('dicodon_distance_matrix.txt', 'w') as f:
    f.write(dicodon_phylip_output)