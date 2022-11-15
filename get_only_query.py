# by Tomás Carrasco
import argparse
from Bio import SeqIO
parser = argparse.ArgumentParser(description= "Script written exclussively to select words from query, and is not a copy of cactus; the header sould have not spaces!!")
parser.add_argument("-word", "--w", help="word to search in the header",type=str)
parser.add_argument("-fasta", "--f", help="fasta",type=str)
arg = parser.parse_args()
word = arg.w
Fasta = open(arg.f)
for rec in SeqIO.parse(Fasta, 'fasta'):
    if word in rec.name:
        print(">"+rec.name)
        print(rec.seq)
