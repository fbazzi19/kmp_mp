Knuth Morris Pratt Algorithm with MultiProcessing

This Python script implements the Knuth Morris Pratt algorithm for pattern searching within a string, utilizing the MultiProcessing library to look for more than one pattern simultaneously.
This program is intended for use specifically with FASTA files containing sequence data and sequence motifs to be searched for within the larger sequence.

Parameters:

python3 kmp_proj.py [seqfile].fasta [patternfile].fasta  
[seqfile].fasta : FASTA file containing the larger sequence that patterns are being searched for in   
[patternfile].fasta : FASTA file containing multiple short sequences to search for in [seqfile].fasta  
