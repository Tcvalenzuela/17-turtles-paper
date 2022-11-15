#python script by Tomás Carrasco

import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description= "filtering the output of orthofinder")
parser.add_argument("-OrthoGroup", "--OG", help="Orthogroup files",type=str)
parser.add_argument("-List", "--L", help="List",type=str)
parser.add_argument("-Output", "--O", help="Output",type=str)

arg = parser.parse_args()
OrthoFinder = open(arg.OG)
List = open(arg.L)
owriter = open(arg.O,"w")
HashList={}
for line in List:
    HashList[line.rstrip("\n")]=0

for line in OrthoFinder:
        if not line.startswith("OG"):
                continue
        else:
                splitter1=line.split("\t")
                OG=splitter1[0]
                splitter2=splitter1[1].split(";")
        for ele in splitter2:
                if ele.startswith("gene"):
                        gene=ele.lstrip("gene=")
                if ele.startswith("protein"):
                        prot=ele
        if OG in HashList:
                owriter.write(str(OG)+"\t"+str(gene)+"\t"+str(prot)+"\n")
