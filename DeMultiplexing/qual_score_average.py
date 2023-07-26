#!/usr/bin/env python
import bioinfo
import argparse

def get_args():
    parser = argparse.ArgumentParser(description="collect variables for kmer retreval")
    parser.add_argument("-f", "--file", help="input file name", type=str, required=True)
    parser.add_argument("-l", "--length", help="input length of sequence", type=int, required=True)
    return parser.parse_args()
args = get_args()

file=args.file #input file
seqLength=args.length #input file
bases=[] #will hold
for spot in range(seqLength):
    bases.append(0)
with open(file, "r") as fastq:
    l = 0 #line number
    for line in fastq:
        l+=1
        if bioinfo.quality_score_line(l):
            line=line.strip()
            i = 0 #reset index
            for char in line:
                #print(char, i, bioinfo.convert_phred(char))
                bases[i] = bioinfo.convert_phred(char)+bases[i] #sums all char in position i
                #print(bases[i])
                i+=1
    for index in range(i):
        bases[index] = bases[index]/(l/4) #total/number of input = mean
print(bases)
