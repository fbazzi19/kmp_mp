import multiprocessing
import sys
import numpy as np
from Bio import SeqIO
import os


def findpattern(text, pname, pattern, outfile, lock):
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
    
    i=0 #index on text
    j=0 #index on pattern
    while (textlen-i) >= (patternlen-j):  
        if text[i]==pattern[j]: #match
            i+=1
            j+=1

        if j==patternlen: #reached the end of the pattern
            lock.acquire()
            #print("I am here w seq " + pname)
            outfile.write("Sequence "+ pname +" found at location "+ str((i-j)+1) +"-"+str((i-j)+patternlen)+"\n")
            print("Seq "+pname+" just wrote to file")
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
    #take the fasta files as arguments when the program is run
    seqfile = sys.argv[1]
    patternsfile = sys.argv[2] #TODO: add a sanity check to return an error if the files arent fasta

    seq = SeqIO.parse(open(seqfile),'fasta')
    for fasta in seq:
        #obtain the sequence as a string from the fasta file
        name, sequence = fasta.id, str(fasta.seq)

    patterns = SeqIO.parse(open(patternsfile),'fasta') 
    f = open("kmpout.txt", "w") #open file to write the location of the patterns to
    proclist=list()#list to hold the processes for each pattern
    for fasta in patterns:
        #obtain each pattern from the fasta file and send it to the find pattern function
        name, pattern = fasta.id, str(fasta.seq)
        #findpattern(sequence, name, pattern, f)
        lock=multiprocessing.Lock()
        px= multiprocessing.Process(target=findpattern, args=(sequence, name, pattern, f, lock))
        proclist.append(px)
        px.start()
        #print("ID of process p1: {}".format(px.pid)) 
    
    for px in proclist:
        px.join()

    f.close() #close output file

    
