#Script by Tomás Carrasco-Valenzuela

import argparse
from Bio import SeqIO
parser = argparse.ArgumentParser(description= "From fasta files print independdant Scaffs")
parser.add_argument("--fa", "-fasta", help="Fasta file to be trimmed")
parser.add_argument("--scaff", "-scaff", help="Scaffolds to output")
arg = parser.parse_args()
Scaffolds=arg.scaff
Hash={}
List={}
openlist=open(arg.scaff)
for line in openlist:
  List[line]=line.rstrip("\n")
inputopen = SeqIO.parse(open(arg.fa),"fasta")
for line in inputopen:
    name=line.id
    seq=line.seq
    splitted=line.id.split(";") #gene;gene=thisthatyouwant;andmoreinfo
    Scaff=splitted[2]
    geneName=Scaff.split("=")[1]
    for ele in List:
      if geneName in ele:
        if name in Hash:
          continue
        else: 
          Hash[name]=(seq)

for rec in Hash:
  print(str(">")+str(rec)+str("\n")+str(Hash[rec]))

#To run:

#python filtering_fasta_from_list.py -fasta GCF_013100865.1_CAS_Tse_1.0translated_cds_long_renamed.fasta -scaff lost_genes_dc.txt > lost_genes_dc_sequences.fasta
