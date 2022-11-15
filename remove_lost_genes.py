# by Tomás Carrasco
import argparse
from Bio import SeqIO
parser = argparse.ArgumentParser(description= "Stolen from myself at From the out of Cactus print in file the chain format")
parser.add_argument("-table", "--t", help="table of genes to delete",type=str)
parser.add_argument("-file", "--f", help="file with genes",type=str)
arg = parser.parse_args()
delnames = open(arg.t)
Genes = open(arg.f)
delete={}
for ele in delnames:
    delete[ele.rstrip("\n")]=1
Keeper={}
NottoKeep={}
for rec in SeqIO.parse(Genes, 'fasta'):

        if rec.name in delete:
            NottoKeep[rec.name]=0
        if rec.name not in delete:   #if rec.name not in delete:
            Keeper[rec.name]=rec.seq
for ele in Keeper:
    print (">"+ele)
    print (Keeper[ele])
