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
    #calculate the length of the text and the pattern
    patternlen = len(pattern)
    textlen = len(text)

    #create an array of zeros to hold lps
    lps = np.zeros((patternlen,), dtype=int)

    #calculate lps array
    lps = computelps(pattern, patternlen, lps)
    print(lps)
    return 0

def computelps(pattern, patternlen, lps):
    """
    Inputs:
    pattern: the pattern to be preprocessed
    patternlen: the length of the pattern
    lps: an array of zeros of size patternlen
    Outputs:
    lps: an array of size pattern len containing the previous location of a prefix within the pattern, if any
    Function:
    Fills an array based on the location of prefixes and reoccuring patterns within the pattern string
    """
    i=1 #idx on the location within pattern
    len=0 #idx on the length of the matching prefix

    while i<patternlen: #increment through the pattern
        if pattern[len]==pattern[i]: #if the location matches a prefix
            len+=1
            lps[i]=len
            i+=1
        else:
            if len>0: #if the length was greather than zero and it doesn't match
                len=lps[len-1]
            else: #if the length was zero and it doesn't match
                lps[i]=0
                i+=1


    return lps


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

