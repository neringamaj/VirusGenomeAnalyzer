<h1>Codon and Dicodon Frequency Analysis in Viral Genomes</h1>

<h2>Repository Description</h2>
<p>
    This repository contains a bioinformatics tool developed to analyze the frequency differences of codons and dicodons in mammalian and bacterial viruses. The sequences, provided in FASTA format, are processed to identify open reading frames, translate nucleotide sequences to amino acids, and calculate the frequency of each codon and dicodon.
</p>

<h2>Features of the Script</h2>
<ul>
    <li><strong>Open Reading Frame Identification:</strong> Identifies all start and stop codon pairs without any intervening stop codons within the provided sequences.</li>
    <li><strong>Sequence Translation:</strong> Translates coding sequences into amino acids, considering codons and dicodons at the protein level.</li>
    <li><strong>Frequency Analysis:</strong> Calculates the frequencies of all possible codons and dicodons within the sequences.</li>
    <li><strong>Distance Matrix Computation:</strong> Computes a distance matrix based on codon and dicodon frequencies to assess the evolutionary relationships between the viruses.</li>
</ul>

<h2>Dependencies</h2>
<ul>
    <li>Python</li>
    <li>BioPython Library</li>
</ul>
