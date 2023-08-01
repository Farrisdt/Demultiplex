#!/usr/bin/env python

import bioinfo
import argparse
import gzip
import matplotlib.pyplot as plt

#take in file name and line length
def get_args():
    parser = argparse.ArgumentParser(description="collect variables for kmer retreval")
    parser.add_argument("-l", "--length", help="hash length (ie. kmer length)", type=int, required=True)
    parser.add_argument("-f", "--file", help="input file name", type=str, required=True)
    parser.add_argument("-o", "--output", help="output file name for histogram", type=str, required=True)
    parser.add_argument("-t", "--filetype", help="Please enter an int between 1 and 4. Options: 1- Read 1, 2- Index 1, 3- Index 2, 4- Read 2", type=int, required=True)
    return parser.parse_args()
args = get_args()

outfilename = args.output
charNum= args.length
currfile=args.file
graphX = [] #for graphing
filetypes = ["Read 1", "Index 1", "Index 2", "Read 2"]

#make list based on given length
scoreList = bioinfo.init_list(charNum)
lineNum = 0 #line counter
#open file
with gzip.open(currfile, "rt") as input:

    for line in input: #go line by line, iterate counter
        lineNum+=1
        line = line.strip()
        if bioinfo.quality_score_line(lineNum):
            for char in range(len(line)): #phred score char and add sum to list
                scoreList[char] += bioinfo.convert_phred(line[char])
#at end divide all sums by number of lines
lineNum=lineNum/4
for base in range(len(scoreList)):
    graphX.append(base+1)
    scoreList[base] = (scoreList[base]/lineNum)

    

#output histogram

plt.bar(graphX, scoreList, width=0.7)
plt.title(f"Distribution of Quality Scores Over Read Location for {filetypes[args.filetype -1]}")
plt.ylabel("Quality Score")
#plt.yscale("log") #makes the y axis on a log scale
plt.xlabel("Base Location on Read")
#plt.xlim([0, 10000]) #max x axis out at 10000
plt.savefig(f"./{outfilename}.png")