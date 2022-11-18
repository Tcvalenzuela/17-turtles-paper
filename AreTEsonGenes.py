import argparse
parser = argparse.ArgumentParser(description= "Identify if TE are in genetic feature and frag the kind of interaction")
parser.add_argument("--TEbed", "-bedfileforTE",  help="bed file with the information of TE")
parser.add_argument("--gff", "-gfffile",  help="Gff file with the genetic features")
arg = parser.parse_args()

TEbed = open(arg.TEbed)
gff = open(arg.gff)
TEStarts={}
TEEnds={}


for TE in TEbed:
    Splitted=TE.rstrip("\n").split("\t")
    Scaff=Splitted[0]
    Start=int(Splitted[1])
    End=int(Splitted[2])
    Kind=str(Splitted[3])
    Kimura=float(Splitted[4])
    dNdS=float(Splitted[5])
    Size=End-Start
    TEStarts[TE]=(Start,Scaff)
    TEEnds[TE]=(End,Scaff)
    TEdNdS[TE]=(dNdS,Kimura)
    TEsize[TE]=(Size,Kind)
Gffexon={}
GffCDS={}
Gffgene={}
GffmRNA={}
Gfftranscript={}
GffcDNA_match={}
for line in gff:
    if line.startswith("#"):
        continue
    Split=line.split("\t")
    gffScaff=Split[0]
    gffclass=Split[2]
    gffStart=int(Split[3])
    gffEnd=int(Split[4])
    if gffclass == "region":
        continue
    elif gffclass == "exon":
        Gffexon[line]=(gffStart,gffEnd,gffScaff)
    elif gffclass == "CDS":
        GffCDS[line]=(gffStart,gffEnd,gffScaff)
    elif gffclass == "mRNA":
        GffmRNA[line]=(gffStart,gffEnd,gffScaff)
    elif gffclass == "gene":
        Gffgene[line]=(gffStart,gffEnd,gffScaff)
    elif gffclass == "transcript":
        Gfftranscript[line]=(gffStart,gffEnd,gffScaff)
    elif gffclass == "cDNA_match":
        GffcDNA_match[line]=(gffStart,gffEnd,gffScaff)
for ele in Gffgene:
    for rec in TEStarts:
        if Gffgene[ele][2]!= TEStarts[rec][1]:
            continue
        elif Gffgene[ele][0] > TEStarts[rec][0] and Gffgene[ele][0] < TEEnds[rec] and  Gffgene[ele][1] < TEEnds[rec] and TEStarts[rec][0] < Gffgene[ele][1]:
            print (rec.rstrip("\n") +"\t"+ "gene" + "_inside_TE"+"\t" +ele + "\t" + TEdNdS[TE][0]+  "\t" +TEdNdS[TE][1]+  "\t" +TEsize[TE][0]+ "\t" +TEsize[TE][0])
        elif Gffgene[ele][0] < TEStarts[rec][0] and Gffgene[ele][0] < TEEnds[rec] and  Gffgene[ele][1] > TEEnds[rec] and TEStarts[rec][0] < Gffgene[ele][1]:
            print (rec.rstrip("\n") + "\t"+"TE_inside_"+ "gene" +"\t"+ ele+ TEdNdS[TE][0]+  "\t" +TEdNdS[TE][1]+  "\t" +TEsize[TE][0]+ "\t" +TEsize[TE][0])
        elif Gffgene[ele][0] <= TEStarts[rec][0] and Gffgene[ele][0] < TEEnds[rec] and  Gffgene[ele][1] <= TEEnds[rec] and TEStarts[rec][0] < Gffgene[ele][1]:
            print (rec.rstrip("\n") +"\t"+ "TE_fromInside2Downstream_"+ "gene" +"\t"+ ele+ TEdNdS[TE][0]+  "\t" +TEdNdS[TE][1]+  "\t" +TEsize[TE][0]+ "\t" +TEsize[TE][0])
        elif Gffgene[ele][0] >= TEStarts[rec][0] and Gffgene[ele][0] < TEEnds[rec] and  Gffgene[ele][1] >= TEEnds[rec] and TEStarts[rec][0] < Gffgene[ele][1]:
            print (rec.rstrip("\n") +"\t"+ "TE_fromInside2upstream_"+"gene" +"\t"+ ele+ TEdNdS[TE][0]+  "\t" +TEdNdS[TE][1]+  "\t" +TEsize[TE][0]+ "\t" +TEsize[TE][0])
        elif Gffgene[ele][0] == TEStarts[rec][0] and Gffgene[ele][0] < TEEnds[rec] and  Gffgene[ele][1] == TEEnds[rec] and TEStarts[rec][0] < Gffgene[ele][1] :
            print (rec.rstrip("\n") +"\t"+ "TE_exactly_"+ "gene"+"\t" + ele+ TEdNdS[TE][0]+  "\t" +TEdNdS[TE][1]+  "\t" +TEsize[TE][0]+ "\t" +TEsize[TE][0])
for ele in GffmRNA:
    for rec in TEStarts:
        if GffmRNA[ele][2]!= TEStarts[rec][1]:
            continue
        elif GffmRNA[ele][0] > TEStarts[rec][0] and GffmRNA[ele][0] < TEEnds[rec] and  GffmRNA[ele][1] < TEEnds[rec] and TEStarts[rec][0] < GffmRNA[ele][1]:
            print (rec.rstrip("\n") +"\t"+ "mRNA" + "_inside_TE"+"\t" +ele+ TEdNdS[TE][0]+  "\t" +TEdNdS[TE][1]+  "\t" +TEsize[TE][0]+ "\t" +TEsize[TE][0])
        elif GffmRNA[ele][0] < TEStarts[rec][0] and GffmRNA[ele][0] < TEEnds[rec] and  GffmRNA[ele][1] > TEEnds[rec] and TEStarts[rec][0] < GffmRNA[ele][1]:
            print (rec.rstrip("\n") + "\t"+"TE_inside_"+ "mRNA" +"\t"+ ele+ TEdNdS[TE][0]+  "\t" +TEdNdS[TE][1]+  "\t" +TEsize[TE][0]+ "\t" +TEsize[TE][0])
        elif GffmRNA[ele][0] <= TEStarts[rec][0] and GffmRNA[ele][0] < TEEnds[rec] and  GffmRNA[ele][1] <= TEEnds[rec] and TEStarts[rec][0] < GffmRNA[ele][1]:
            print (rec.rstrip("\n") +"\t"+ "TE_fromInside2Downstream_"+ "mRNA" +"\t"+ ele+ TEdNdS[TE][0]+  "\t" +TEdNdS[TE][1]+  "\t" +TEsize[TE][0]+ "\t" +TEsize[TE][0])
        elif GffmRNA[ele][0] >= TEStarts[rec][0] and GffmRNA[ele][0] < TEEnds[rec] and  GffmRNA[ele][1] >= TEEnds[rec] and TEStarts[rec][0] < GffmRNA[ele][1]:
            print (rec.rstrip("\n") +"\t"+ "TE_fromInside2upstream_"+"mRNA" +"\t"+ ele+ TEdNdS[TE][0]+  "\t" +TEdNdS[TE][1]+  "\t" +TEsize[TE][0]+ "\t" +TEsize[TE][0])
        elif GffmRNA[ele][0] == TEStarts[rec][0] and GffmRNA[ele][0] < TEEnds[rec] and  GffmRNA[ele][1] == TEEnds[rec] and TEStarts[rec][0] < GffmRNA[ele][1] :
            print (rec.rstrip("\n") +"\t"+ "TE_exactly_"+ "mRNA"+"\t" + ele+ TEdNdS[TE][0]+  "\t" +TEdNdS[TE][1]+  "\t" +TEsize[TE][0]+ "\t" +TEsize[TE][0])
for ele in Gffexon:
    for rec in TEStarts:
        if Gffexon[ele][2]!= TEStarts[rec][1]:
            continue
        elif Gffexon[ele][0] > TEStarts[rec][0] and Gffexon[ele][0] < TEEnds[rec] and  Gffexon[ele][1] < TEEnds[rec] and TEStarts[rec][0] < Gffexon[ele][1]:
            print (rec.rstrip("\n") +"\t"+ "exon" + "_inside_TE"+"\t" +ele+ TEdNdS[TE][0]+  "\t" +TEdNdS[TE][1]+  "\t" +TEsize[TE][0]+ "\t" +TEsize[TE][0])
            print (rec.rstrip("\n") + "\t"+"TE_inside_"+ "exon" +"\t"+ ele+ TEdNdS[TE][0]+  "\t" +TEdNdS[TE][1]+  "\t" +TEsize[TE][0]+ "\t" +TEsize[TE][0])
        elif Gffexon[ele][0] <= TEStarts[rec][0] and Gffexon[ele][0] < TEEnds[rec] and  Gffexon[ele][1] <= TEEnds[rec] and TEStarts[rec][0] < Gffexon[ele][1]:
            print (rec.rstrip("\n") +"\t"+ "TE_fromInside2Downstream_"+ "exon" +"\t"+ ele+ TEdNdS[TE][0]+  "\t" +TEdNdS[TE][1]+  "\t" +TEsize[TE][0]+ "\t" +TEsize[TE][0])
        elif Gffexon[ele][0] >= TEStarts[rec][0] and Gffexon[ele][0] < TEEnds[rec] and  Gffexon[ele][1] >= TEEnds[rec] and TEStarts[rec][0] < Gffexon[ele][1]:
            print (rec.rstrip("\n") +"\t"+ "TE_fromInside2upstream_"+"exon" +"\t"+ ele+ TEdNdS[TE][0]+  "\t" +TEdNdS[TE][1]+  "\t" +TEsize[TE][0]+ "\t" +TEsize[TE][0])
        elif Gffexon[ele][0] == TEStarts[rec][0] and Gffexon[ele][0] < TEEnds[rec] and  Gffexon[ele][1] == TEEnds[rec] and TEStarts[rec][0] < Gffexon[ele][1] :
            print (rec.rstrip("\n") +"\t"+ "TE_exactly_"+ "exon"+"\t" + ele+ TEdNdS[TE][0]+  "\t" +TEdNdS[TE][1]+  "\t" +TEsize[TE][0]+ "\t" +TEsize[TE][0])
for ele in GffCDS:
    for rec in TEStarts:
        if GffCDS[ele][2]!= TEStarts[rec][1]:
            continue
        elif GffCDS[ele][0] > TEStarts[rec][0] and GffCDS[ele][0] < TEEnds[rec] and  GffCDS[ele][1] < TEEnds[rec] and TEStarts[rec][0] < GffCDS[ele][1]:
            print (rec.rstrip("\n") +"\t"+ "CDS" + "_inside_TE"+"\t" +ele+ TEdNdS[TE][0]+  "\t" +TEdNdS[TE][1]+  "\t" +TEsize[TE][0]+ "\t" +TEsize[TE][0])
        elif GffCDS[ele][0] < TEStarts[rec][0] and GffCDS[ele][0] < TEEnds[rec] and  GffCDS[ele][1] > TEEnds[rec] and TEStarts[rec][0] < GffCDS[ele][1]:
            print (rec.rstrip("\n") + "\t"+"TE_inside_"+ "CDS" +"\t"+ ele+ TEdNdS[TE][0]+  "\t" +TEdNdS[TE][1]+  "\t" +TEsize[TE][0]+ "\t" +TEsize[TE][0])
        elif GffCDS[ele][0] <= TEStarts[rec][0] and GffCDS[ele][0] < TEEnds[rec] and  GffCDS[ele][1] <= TEEnds[rec] and TEStarts[rec][0] < GffCDS[ele][1]:
            print (rec.rstrip("\n") +"\t"+ "TE_fromInside2Downstream_"+ "CDS" +"\t"+ ele+ TEdNdS[TE][0]+  "\t" +TEdNdS[TE][1]+  "\t" +TEsize[TE][0]+ "\t" +TEsize[TE][0])
        elif GffCDS[ele][0] >= TEStarts[rec][0] and GffCDS[ele][0] < TEEnds[rec] and  GffCDS[ele][1] >= TEEnds[rec] and TEStarts[rec][0] < GffCDS[ele][1]:
            print (rec.rstrip("\n") +"\t"+ "TE_fromInside2upstream_"+"CDS" +"\t"+ ele+ TEdNdS[TE][0]+  "\t" +TEdNdS[TE][1]+  "\t" +TEsize[TE][0]+ "\t" +TEsize[TE][0])
        elif GffCDS[ele][0] == TEStarts[rec][0] and GffCDS[ele][0] < TEEnds[rec] and  GffCDS[ele][1] == TEEnds[rec] and TEStarts[rec][0] < GffCDS[ele][1] :
            print (rec.rstrip("\n") +"\t"+ "TE_exactly_"+ "CDS"+"\t" + ele+ TEdNdS[TE][0]+  "\t" +TEdNdS[TE][1]+  "\t" +TEsize[TE][0]+ "\t" +TEsize[TE][0])
for ele in Gfftranscript:
    for rec in TEStarts:
        if Gfftranscript[ele][2]!= TEStarts[rec][1]:
            continue
        elif Gfftranscript[ele][0] > TEStarts[rec][0] and Gfftranscript[ele][0] < TEEnds[rec] and  Gfftranscript[ele][1] < TEEnds[rec] and TEStarts[rec][0] < Gfftranscript[ele][1]:
            print (rec.rstrip("\n") +"\t"+ "transcript" + "_inside_TE"+"\t" +ele+ TEdNdS[TE][0]+  "\t" +TEdNdS[TE][1]+  "\t" +TEsize[TE][0]+ "\t" +TEsize[TE][0])
        elif Gfftranscript[ele][0] < TEStarts[rec][0] and Gfftranscript[ele][0] < TEEnds[rec] and  Gfftranscript[ele][1] > TEEnds[rec] and TEStarts[rec][0] < Gfftranscript[ele][1]:
            print (rec.rstrip("\n") + "\t"+"TE_inside_"+ "transcript" +"\t"+ ele+ TEdNdS[TE][0]+  "\t" +TEdNdS[TE][1]+  "\t" +TEsize[TE][0]+ "\t" +TEsize[TE][0])
        elif Gfftranscript[ele][0] <= TEStarts[rec][0] and Gfftranscript[ele][0] < TEEnds[rec] and  Gfftranscript[ele][1] <= TEEnds[rec] and TEStarts[rec][0] < Gfftranscript[ele][1]:
            print (rec.rstrip("\n") +"\t"+ "TE_fromInside2Downstream_"+ "transcript" +"\t"+ ele+ TEdNdS[TE][0]+  "\t" +TEdNdS[TE][1]+  "\t" +TEsize[TE][0]+ "\t" +TEsize[TE][0])
        elif Gfftranscript[ele][0] >= TEStarts[rec][0] and Gfftranscript[ele][0] < TEEnds[rec] and  Gfftranscript[ele][1] >= TEEnds[rec] and TEStarts[rec][0] < Gfftranscript[ele][1]:
            print (rec.rstrip("\n") +"\t"+ "TE_fromInside2upstream_"+"transcript" +"\t"+ ele+ TEdNdS[TE][0]+  "\t" +TEdNdS[TE][1]+  "\t" +TEsize[TE][0]+ "\t" +TEsize[TE][0])
        elif Gfftranscript[ele][0] == TEStarts[rec][0] and Gfftranscript[ele][0] < TEEnds[rec] and  Gfftranscript[ele][1] == TEEnds[rec] and TEStarts[rec][0] < Gfftranscript[ele][1] :
            print (rec.rstrip("\n") +"\t"+ "TE_exactly_"+ "transcript"+"\t" + ele+ TEdNdS[TE][0]+  "\t" +TEdNdS[TE][1]+  "\t" +TEsize[TE][0]+ "\t" +TEsize[TE][0])
for ele in GffcDNA_match:
    for rec in TEStarts:
        if GffcDNA_match[ele][2]!= TEStarts[rec][1]:
            continue
        elif GffcDNA_match[ele][0] > TEStarts[rec][0] and GffcDNA_match[ele][0] < TEEnds[rec] and  GffcDNA_match[ele][1] < TEEnds[rec] and TEStarts[rec][0] < GffcDNA_match[ele][1]:
            print (rec.rstrip("\n") +"\t"+ "cDNA_match" + "_inside_TE"+"\t" +ele+ TEdNdS[TE][0]+  "\t" +TEdNdS[TE][1]+  "\t" +TEsize[TE][0]+ "\t" +TEsize[TE][0])
        elif GffcDNA_match[ele][0] < TEStarts[rec][0] and GffcDNA_match[ele][0] < TEEnds[rec] and  GffcDNA_match[ele][1] > TEEnds[rec] and TEStarts[rec][0] < GffcDNA_match[ele][1]:
            print (rec.rstrip("\n") + "\t"+"TE_inside_"+ "cDNA_match" +"\t"+ ele+ TEdNdS[TE][0]+  "\t" +TEdNdS[TE][1]+  "\t" +TEsize[TE][0]+ "\t" +TEsize[TE][0])
        elif GffcDNA_match[ele][0] <= TEStarts[rec][0] and GffcDNA_match[ele][0] < TEEnds[rec] and  GffcDNA_match[ele][1] <= TEEnds[rec] and TEStarts[rec][0] < GffcDNA_match[ele][1]:
            print (rec.rstrip("\n") +"\t"+ "TE_fromInside2Downstream_"+ "cDNA_match" +"\t"+ ele+ TEdNdS[TE][0]+  "\t" +TEdNdS[TE][1]+  "\t" +TEsize[TE][0]+ "\t" +TEsize[TE][0])
        elif GffcDNA_match[ele][0] >= TEStarts[rec][0] and GffcDNA_match[ele][0] < TEEnds[rec] and  GffcDNA_match[ele][1] >= TEEnds[rec] and TEStarts[rec][0] < GffcDNA_match[ele][1]:
            print (rec.rstrip("\n") +"\t"+ "TE_fromInside2upstream_"+"cDNA_match" +"\t"+ ele+ TEdNdS[TE][0]+  "\t" +TEdNdS[TE][1]+  "\t" +TEsize[TE][0]+ "\t" +TEsize[TE][0])
        elif GffcDNA_match[ele][0] == TEStarts[rec][0] and GffcDNA_match[ele][0] < TEEnds[rec] and  GffcDNA_match[ele][1] == TEEnds[rec] and TEStarts[rec][0] < GffcDNA_match[ele][1] :
            print (rec.rstrip("\n") +"\t"+ "TE_exactly_"+ "cDNA_match"+"\t" + ele+ TEdNdS[TE][0]+  "\t" +TEdNdS[TE][1]+  "\t" +TEsize[TE][0]+ "\t" +TEsize[TE][0])

