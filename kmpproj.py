import multiprocessing
import sys
import numpy as np
from Bio import SeqIO

def findpattern(text, pname, pattern):
    """
    Inputs:
    text : string pattern is being searched for in
    pname : the ID of the pattern found
    pattern : string being searched for
    Outputs:
    none
    Function
    Searches for pattern in text and writes the pattern ID and the location of the pattern to a text file
    """
    return 0

def computelps(pattern, patternlen, lps):
    return 0


if __name__=="__main__":
    #take the fasta files as arguments when the program is run
    seqfile = sys.argv[1]
    patternsfile = sys.argv[2] #TODO: add a sanity check to return an error if the files arent fasta

    seq = SeqIO.parse(open(seqfile),'fasta')
    for fasta in seq:
        #obtain the sequence as a string from the fasta file
        name, sequence = fasta.id, str(fasta.seq)

    patterns = SeqIO.parse(open(patternsfile),'fasta') #TO CONSIDER: Should the writing to output file happen here or within function
    for fasta in patterns:
        #obtain each pattern from the fasta file and send it to the find pattern function
        name, pattern = fasta.id, str(fasta.seq)
        findpattern(sequence, name, pattern)

