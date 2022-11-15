# by Tomás Carrasco
from Bio import SeqIO
import sys
input = sys.argv[1] # python command to recognize the first argument as the fasta file I want to open

hash={} # open a table to store GeneID
fasta_sequences = SeqIO.parse(open(input),"fasta") #read every sequence on fasta file
for rec in fasta_sequences: 
        stripper=rec.id.split("_") #split the header of the sequence and identify the number of the GeneID that apears after ':' in the headers
        if stripper[0] in hash: #starts the looṕ, everytime it reads a GeneID it will store it on hash
                if len(rec.seq) > len(hash[stripper[0]].seq): #look at sequence by sequence and if the sequence of an GeneID were longer
                        hash[stripper[0]]=rec # and if the GeneID is already in the hash, it will store the new sequence
        else:
                hash[stripper[0]]=rec #store in hash every sequence which GeneID appears for the first and only time
for rec in hash:
        print(">"+hash[rec].id) #print everything on record for the hash table
        print(hash[rec].seq)
