# Knuth Morris Pratt Algorithm with MultiProcessing

This Python script implements the Knuth Morris Pratt algorithm for pattern searching within a string, utilizing the MultiProcessing library to look for more than one pattern simultaneously.
This program is intended for use specifically with FASTA files containing sequence data and sequence motifs to be searched for within the larger sequence.

## Parameters:

`python3 kmp_proj.py [seqfile].fasta [patternfile].fasta`  
[seqfile].fasta: FASTA file containing the larger sequence that patterns are being searched for in   
[patternfile].fasta: FASTA file containing multiple short sequences to search for in [seqfile].fasta  

### Constraints  
There must be two parameters provided.  
Files must be FASTA files.  
The FASTA file containing the larger sequence must only have one sequence.  
The FASTA file containing the patterns must be specified after the larger sequence file.  

## Output:  
`kmpout.txt`  
A text file containing several lines in the format:  
> Sequence `[sequence name]` found at location `[start index]`-`[end index]`

Specifying every location where a pattern was found.  

## Test Files


## References  
For more information about the Knuth Morris Pratt algorithm, see [the Wikipedia page](https://en.wikipedia.org/wiki/Knuth%E2%80%93Morris%E2%80%93Pratt_algorithm#). Much thanks to [this YouTube video](https://www.youtube.com/watch?v=V5-7GzOfADQ) also. FASTA files containing sequence data from chromosome 7 of the human genome along with the genome of SAR-COV-2 gathered from the [National Center for Biotechnology Information](https://www.ncbi.nlm.nih.gov/).
