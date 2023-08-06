#!/usr/bin/env python
import argparse
import bioinfo
import numpy as np
import gzip
from matplotlib import pyplot as plt

def get_args():
    parser = argparse.ArgumentParser(description="collect variables for kmer retreval")
    parser.add_argument("-r1", "--read1", help="read 1 file name", type=str, required=True)
    parser.add_argument("-r2", "--read2", help="read 2 file name", type=str, required=True)
    parser.add_argument("-i1", "--index1", help="index 1 file name", type=str, required=True)
    parser.add_argument("-i2", "--index2", help="index 2 file name", type=str, required=True)        
    parser.add_argument("-i", "--indexL", help="index length", default=8, required=False)
    parser.add_argument("-q", "--qualcut", help="quality score cutoff", default=30, required=False)
    parser.add_argument("-l", "--lineNumber", help="number of lines in the file", default=1452986940, required=False)
    return parser.parse_args()
args = get_args()

#Variables
r1file = args.read1
r2file = args.read2
i1file = args.index1
i2file = args.index2
qualcutoff = args.qualcut #quality cutoff score, use argparse
indexLength = args.indexL #argparse
hasN = False # to check for 'N' in index
lineNum= int(args.lineNumber) #numbner of lines in file, found in data exploration.
record1 = "" # for writting out
record2 = ""
r1seq = "" # for curr index
r2seq = ""
indexs = {"GTAGCGTA": "B1", "CGATCGAT": "A5", "GATCAAGG": "C1", "AACAGCGA": "B9", "TAGCCATG": "C9", "CGGTAATC": "C3", "CTCTGGAT": "B3", "TACCGGAT": "C4", "CTAGCTCA": "A11", "CACTTCAC": "C7", "GCTACTCT": "B2", "ACGATCAG": "A1", "TATGGCAC": "B7", "TGTTCCGT": "A3", "GTCCTAAG": "B4", "TCGACAAG": "A12", "TCTTCGAC": "C10", "ATCATGCG": "A2", "ATCGTGGT": "C2", "TCGAGAGT": "A10", "TCGGATTC": "B8", "GATCTTGC": "A7", "AGAGTCCA": "B10", "AGGATAGC": "A8"}
heatmap = {}
for keys in indexs:
    indexs[keys] = [0, 0, 0] #times properly matched, times mismatched in read 1, times mismatched in read 2
    for key2 in indexs:
        input = f"{keys}:{key2}"
        heatmap[input] = 0 #Holds mismatch encounter rates, keys is on r1 key2 is on r2.
mismatch = 0 # mismatch total
indexs["Unknown"]= [0, 0, 0] # times both indexs unknown, times only read 1 unknown, times only read 2 unknown
files = {} #to hold open file names

#open all 52 output files
#maybe add in something to prevent overwrite
r1mismatched = open("./output/r1_Mismatched.fastq", "w")
r2mismatched = open("./output/r2_Mismatched.fastq", "w")
for key in indexs:
    r1name=f"./output/r1_{key}.fastq"
    r2name=f"./output/r2_{key}.fastq"
    files[key] = (open(r1name, "w"), open(r2name, "w"))

with gzip.open(r1file, "rt") as r1, gzip.open(i1file, "rt") as i1, gzip.open(i2file, "rt") as i2, gzip.open(r2file, "rt") as r2:
#with open(r1file, "r") as r1, open(i1file, "r") as i1, open(i2file, "r") as i2, open(r2file, "r") as r2: #for testing
    i=0
    while i<lineNum:
        i+=1
        a = (r1.readline()).strip()
        b = (i1.readline()).strip()
        craw = (i2.readline()).strip()
        d = (r2.readline()).strip()
        if bioinfo.header_line(i): #make bioinfo script to check line number based on i, return true if header
            h1 = a
            h2 = d
        if bioinfo.DNA_line(i): #same as above, may already have coded
            c = bioinfo.reverse_complement(craw) #make for bioinfo. Turn A>T etc and reverse order, return new string . N>N. toupper everything. error not DNA if not ACTGN.
            r1seq = b
            r2seq = c
            h1+=(f" {b}-{c}")
            record1=f"{a}\n+\n"
            h2+=(f" {b}-{c}")
            record2=f"{d}\n+\n"
        if bioinfo.quality_score_line(i): #already in bioinfo
            record1+=f"{a}\n"
            record2+=f"{d}\n"
            hasN = False
            for char in b:
                if char=='N':
                    hasN=True
                    files["Unknown"][0].write(f"{h1}\n{record1}")
                    indexs["Unknown"][1] = indexs["Unknown"][1] + 1
                    break
            for char in c:
                if char=='N':
                    if(not hasN):
                        files["Unknown"][1].write(f"{h2}\n{record2}")
                        indexs["Unknown"][0] = indexs["Unknown"][0] + 1
                    indexs["Unknown"][2] = indexs["Unknown"][2] + 1
                    break
            if bioinfo.qual_score(b)<qualcutoff or bioinfo.qual_score(c)<qualcutoff:
                files["Unknown"][0].write(f"{h1}\n{record1}")
                files["Unknown"][1].write(f"{h2}\n{record2}")               
            elif r1seq in indexs and r2seq in indexs:
                if r1seq==r2seq:
                    files[r1seq][0].write(f"{h1}\n{record1}")
                    files[r2seq][1].write(f"{h2}\n{record2}")
                    indexs[r1seq][0] = indexs[r1seq][0] + 1
                else:
                    indexs[r1seq][1] = indexs[r1seq][1] + 1
                    indexs[r2seq][2] = indexs[r2seq][2] + 1
                    input = f"{r1seq}:{r2seq}"
                    heatmap[input]+=1 
                    mismatch = mismatch + 1
                    r1mismatched.write(f"{h1}\n{record1}")
                    r2mismatched.write(f"{h2}\n{record2}")
            else:
                files["Unknown"][0].write(f"{h1}\n{record1}")
                files["Unknown"][1].write(f"{h2}\n{record2}")
                indexs["Unknown"][0] = indexs["Unknown"][0] + 1
                if not (r1seq in indexs):
                    indexs["Unknown"][1] = indexs["Unknown"][1] + 1
                if not (r2seq in indexs):
                    indexs["Unknown"][2] = indexs["Unknown"][2] + 1
#close all 52 output files
r1mismatched.close()
r2mismatched.close()
for key in indexs:
    files[key][0].close()
    files[key][1].close()

output = open("./output/report.txt", "w")

matchTotal=0
for key in indexs:
    matchTotal+=indexs[key][0]
output.write(f"Total Matched: {matchTotal}\n")
output.write(f'Total Mismatched: {mismatch}\n')
output.write("Index Matched Mismatched:Read 1  Mismatched:Read 2\n")
for key in sorted(indexs): #, key=indexsitems()
    output.write(f"{key}\t{indexs[key][0]}\t{indexs[key][1]}\t{indexs[key][2]}\n")

#graph
graph = np.zeros((24,24))
x= -1
y= -1 
xAxis=[]
yAxis=[]
for key1 in sorted(indexs):
    x+=1
    if not (key1=='Unknown'):
        xAxis.append(key1)
        for key2 in sorted(indexs):
            if not (key2=='Unknown'):
                y+=1
                graph[x,y] = heatmap[f"{key1}:{key2}"]
        yAxis.append(key2)
    y=-1
for item in xAxis:
    output.write(f"\t{item}") #header
output.write("\n")
curr=0
for item in xAxis:
    output.write(f"{item}\t") #header
    for i in range(len(xAxis)):
        amount=graph[i,curr]
        percent=(amount/mismatch)*100
        output.write(f"{amount}({percent}%)\t") #data
    curr+=1
    output.write("\n")

output.close

'''Graphing Graveyard'''
#plt.table(graph)
#plt.imshow(graph, cmpa='hot', interpolation='nearest')

#dataTable = plt.table(cellText=graph, rowLabels=xAxis, colLabels=yAxis, loc='top')


'''
fig, axs =plt.subplots(24,24)
clust_data = np.random.random((24,24))
#collabel=("col 1", "col 2", "col 3")
#axs[0].axis('tight')
#axs[0].axis('off')
the_table = axs[0].table(cellText=clust_data,colLabels=xAxis,loc='center')
#axs[1].plot(clust_data[:,0],clust_data[:,1])
'''

#plt.savefig('TableandHeatmap.png')