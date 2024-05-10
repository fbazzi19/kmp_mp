import multiprocessing
import sys
import numpy as np
from Bio import SeqIO
import os


def findpattern(text, pname, pattern, mpqueue, lock):
    """
    Inputs:
    text : string pattern is being searched for in
    pname : the ID of the pattern found
    pattern : string being searched for
    Outputs:
    none
    Function:
    Searches for pattern in text and writes the pattern ID and the location of the pattern to a multiprocessing queue.
    The location is indicative of the index within text, rather than actual chromosome location. Locations are indexed 
    starting at 1. For more information on the Knuth Morris Pratt algorith used to search for the pattern, see
    https://en.wikipedia.org/wiki/Knuth%E2%80%93Morris%E2%80%93Pratt_algorithm#
    """
    #calculate the length of the text and the pattern
    patternlen = len(pattern)
    textlen = len(text)

    #create an array of zeros to hold lps
    lps = np.zeros((patternlen,), dtype=int)

    #calculate lps array
    lps = computelps(pattern, patternlen, lps)
    
    i=0 #index on text
    j=0 #index on pattern
    while (textlen-i) >= (patternlen-j):  
        if text[i]==pattern[j]: #match
            i+=1
            j+=1

        if j==patternlen: #reached the end of the pattern
            lock.acquire()
            mpqueue.put("> "+ pname +": "+ str((i-j)+1) +"-"+str((i-j)+patternlen)+"\n") #put the results in a multiprocessing queue
            lock.release()
            j=lps[j-1]

        elif i<textlen and text[i]!=pattern[j]: #not a match
            if j>0:#hasn't backtracked to the start of the pattern
                j=lps[j-1]
            else: #has reached the start of the pattern
                i+=1

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
    #check if there are not enough or too many parameters and return an error if so
    if len(sys.argv)!=3:
        sys.exit("User must specify exactly 2 parameters.")

    #take the fasta files as arguments when the program is run
    seqfile = sys.argv[1]
    patternsfile = sys.argv[2]

    #check if both paramaters are fasta files and return an error if not
    if seqfile.split(".")[-1] != "fasta" or patternsfile.split(".")[-1] != "fasta":
        sys.exit("Both files must be FASTA files.")

    seq = SeqIO.parse(open(seqfile),'fasta')
    numseqs = 0 #count the number of sequences in the text file
    for fasta in seq:
        #check if there is more than one sequence in the text file and return an error if so
        if numseqs > 0:
            sys.exit("The first FASTA file can only contain one sequence.")
        #obtain the sequence as a string from the fasta file
        name, sequence = fasta.id, str(fasta.seq)
        numseqs+=1

    patterns = SeqIO.parse(open(patternsfile),'fasta') 
    
    #queue to put results from all the processes into
    resultsqueue=multiprocessing.Queue()
    #lock to avoid adding to queue conflicts
    lock=multiprocessing.Lock()
    f = open("kmpout.txt", "w") #open file to write the location of the patterns to
    #process for each pattern
    px= [multiprocessing.Process(target=findpattern, args=(sequence, fasta.id, str(fasta.seq), resultsqueue, lock))for fasta in patterns]
    
    #start processes
    for p in px:
        p.start()

    #join processes
    for p in px:
        p.join()

    #write queue to output file
    for i in range(resultsqueue.qsize()):
        f.write(resultsqueue.get())

    f.close() #close output file

    
