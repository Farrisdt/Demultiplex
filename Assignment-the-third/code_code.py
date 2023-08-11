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
    parser.add_argument("-i", "--indexL", help="index length", default=8, type=int, required=False)
    parser.add_argument("-q", "--qualcut", help="quality score cutoff", default=30, type=int, required=False)
    parser.add_argument("-l", "--lineNumber", help="number of lines in the file", default=1452986940, type=int, required=False)
    return parser.parse_args()
args = get_args()

#Variables
r1file = args.read1
r2file = args.read2
i1file = args.index1
i2file = args.index2
qualcutoff = args.qualcut
indexLength = args.indexL
hasN = False # set to true if an index has an N, to prevent poor quality data being written to the unknown file multiple times
lineNum= int(args.lineNumber) #found in data exploration.
record1 = "" # for writting out
record2 = ""
r1seq = "" # for curr index
r2seq = ""
#dict of all expected barcodes along with their well plate location. Well plate location not currently used.
indexs = {"GTAGCGTA": "B1", "CGATCGAT": "A5", "GATCAAGG": "C1", "AACAGCGA": "B9", "TAGCCATG": "C9", "CGGTAATC": "C3", "CTCTGGAT": "B3", "TACCGGAT": "C4", "CTAGCTCA": "A11", "CACTTCAC": "C7", "GCTACTCT": "B2", "ACGATCAG": "A1", "TATGGCAC": "B7", "TGTTCCGT": "A3", "GTCCTAAG": "B4", "TCGACAAG": "A12", "TCTTCGAC": "C10", "ATCATGCG": "A2", "ATCGTGGT": "C2", "TCGAGAGT": "A10", "TCGGATTC": "B8", "GATCTTGC": "A7", "AGAGTCCA": "B10", "AGGATAGC": "A8"}
heatmap = {} #for graphing
for keys in indexs:
    indexs[keys] = [0, 0, 0] #times properly matched, times mismatched in read 1, times mismatched in read 2
    for key2 in indexs:
        input = f"{keys}:{key2}"
        heatmap[input] = 0 #Holds mismatch encounter rates, keys is on r1 key2 is on r2.
mismatch = 0 # mismatch total
indexs["Unknown"]= [0, 0, 0] # the number of records in unknown, times read 1 was unknown, times read 2 was unknown. read 1 and read 2 will have overlap
files = {} #to hold open file names

#open all 52 output files
#currently will overwrite any data in output if run again. Rename or move data that should be kept.
outputname1=f"./output_{qualcutoff}/r1_Mismatched.fastq"
outputname2=f"./output_{qualcutoff}/r2_Mismatched.fastq"
r1mismatched = open(outputname1, "w")
r2mismatched = open(outputname2, "w")
for key in indexs:
    r1name=f"./output_{qualcutoff}/r1_{key}.fastq"
    r2name=f"./output_{qualcutoff}/r2_{key}.fastq"
    files[key] = (open(r1name, "w"), open(r2name, "w"))

with gzip.open(r1file, "rt") as r1, gzip.open(i1file, "rt") as i1, gzip.open(i2file, "rt") as i2, gzip.open(r2file, "rt") as r2:
#with open(r1file, "r") as r1, open(i1file, "r") as i1, open(i2file, "r") as i2, open(r2file, "r") as r2: #for testing with unzipped files
    i=0
    while i<lineNum:
        i+=1
        a = (r1.readline()).strip() #read 1
        b = (i1.readline()).strip() #index 1
        craw = (i2.readline()).strip() #index 2, needs to be reverse complimented
        d = (r2.readline()).strip() #read 2
        if bioinfo.header_line(i): #returns true if header
            h1 = a #read 1 header
            h2 = d #read 2 header
        if bioinfo.DNA_line(i): #for DNA sequence line
            c = bioinfo.reverse_complement(craw) #reverse compliments index 3
            r1seq = b #index 1
            r2seq = c #index 2
            h1+=(f" {b}-{c}") #append indexs to read 1 header
            record1=f"{a}\n+\n" #create record for read 1 output, DNA sequecne and '+' line
            h2+=(f" {b}-{c}") #append indexs to read 2 header
            record2=f"{d}\n+\n" #create record for read 2 output, DNA sequecne and '+' line
        if bioinfo.quality_score_line(i): # if at end of record
            record1+=f"{a}\n" #add quality score to record for read 1
            record2+=f"{d}\n" #add quality score to record for read 2
            hasN = False #reset for loop
            for char in b: #for index 1
                if char=='N': #if any nts were unreadable
                    hasN=True
                    files["Unknown"][0].write(f"{h1}\n{record1}") #write record to read 1 file
                    files["Unknown"][1].write(f"{h2}\n{record2}") #write record to read 2 file
                    indexs["Unknown"][1] = indexs["Unknown"][1] + 1 #count unknown in read 1
                    indexs["Unknown"][0] = indexs["Unknown"][0] + 1 # add to total unknown
                    break #only need 1 N, no need to continue, prevents recording twice
            for char in c: #for index 2
                if char=='N':
                    if(not hasN): # to prevent counting/writing twice
                        hasN=True 
                        files["Unknown"][0].write(f"{h1}\n{record1}") #write record to read 1 file
                        files["Unknown"][1].write(f"{h2}\n{record2}") #write record to read 2 file
                        indexs["Unknown"][0] = indexs["Unknown"][0] + 1 # add to total unknown
                    indexs["Unknown"][2] = indexs["Unknown"][2] + 1 # count unknown in read 2
                    break #only need 1 N, no need to continue, prevents recording twice
            if (bioinfo.qual_score(b)<qualcutoff or bioinfo.qual_score(c)<qualcutoff) and not hasN: # if either index is below the cutoff and the record has not already been recorded
                files["Unknown"][0].write(f"{h1}\n{record1}") #write to read 1
                files["Unknown"][1].write(f"{h2}\n{record2}") #write to read 2
                indexs["Unknown"][0] = indexs["Unknown"][0] + 1 #count total unknown
                if bioinfo.qual_score(b)<qualcutoff:
                    indexs["Unknown"][1] = indexs["Unknown"][1] + 1 #count read 1 unknown
                if bioinfo.qual_score(c)<qualcutoff: 
                    indexs["Unknown"][2] = indexs["Unknown"][2] + 1 #count read 2 unknown     
            elif r1seq in indexs and r2seq in indexs: # if a known index
                if r1seq==r2seq: #if indexs match
                    files[r1seq][0].write(f"{h1}\n{record1}") #write to read 1
                    files[r2seq][1].write(f"{h2}\n{record2}") #write to read 2
                    indexs[r1seq][0] = indexs[r1seq][0] + 1 # count matched index, only need to do once since identical
                else:
                    indexs[r1seq][1] = indexs[r1seq][1] + 1 #count mismatch for current index on read 1
                    indexs[r2seq][2] = indexs[r2seq][2] + 1 #count mismatch for current index on read 2
                    input = f"{r1seq}:{r2seq}" #key for graphing dict
                    heatmap[input]+=1 
                    mismatch = mismatch + 1 #count total mismatch
                    r1mismatched.write(f"{h1}\n{record1}") #write to read 1
                    r2mismatched.write(f"{h2}\n{record2}") #write to read 2
            else:
                if not hasN: #checks to see if already writen from having an N. In line with quality check so no need to confirm that
                    files["Unknown"][0].write(f"{h1}\n{record1}") #write to read 1
                    files["Unknown"][1].write(f"{h2}\n{record2}") #write to read 2
                    indexs["Unknown"][0] = indexs["Unknown"][0] + 1 #count total unknown
                if not (r1seq in indexs):
                    indexs["Unknown"][1] = indexs["Unknown"][1] + 1 #count read 1 unknown
                if not (r2seq in indexs):
                    indexs["Unknown"][2] = indexs["Unknown"][2] + 1 #count read 2 unknown

#close all 52 output files
r1mismatched.close()
r2mismatched.close()
for key in indexs:
    files[key][0].close()
    files[key][1].close()
reportname=f"./output_{qualcutoff}/report.txt"
output = open(reportname, "w")

matchTotal=0
for key in indexs:
    matchTotal+=indexs[key][0] #counting total for all keys
matchTotal=matchTotal-indexs["Unknown"][0] #removing unknown from total
#data output
output.write(f"Quality cut-off: {qualcutoff}\n")
output.write(f"Total Matched: {matchTotal}\n")
output.write(f'Total Mismatched: {mismatch}\n')
output.write(f'Total Unknown: {indexs["Unknown"][0]}\n')
output.write(f'Unknown from read 1: {indexs["Unknown"][1]}\tUnknown from read 2: {indexs["Unknown"][2]}\n\n')
output.write("Index\tMatched\tMismatched:Read 1\tMismatched:Read 2\n") #header
for key in sorted(indexs):
    if not(key=="Unknown"): #unknown reported above
        output.write(f"{key}\t{indexs[key][0]}\t{indexs[key][1]}\t{indexs[key][2]}\n")
output.write("\n") #for readability

#graph
graph = np.zeros((24,24)) #emtpy numpy array
x= -1 #line counters
y= -1 
Axis=[] #holds ordered axis, same for X and Y
for key1 in sorted(indexs):
    x+=1
    if not (key1=='Unknown'): #ignore unknown
        Axis.append(key1) #keeps the same order as the data
        for key2 in sorted(indexs):
            if not (key2=='Unknown'): #ignore unknown
                y+=1
                graph[x,y] = heatmap[f"{key1}:{key2}"] #creates a location sensitive pairing in numpy
    y=-1
for item in Axis:
    output.write(f"\t{item}") #header
output.write("\n")
curr=0
for item in Axis:
    output.write(f"{item}\t") #y axis, one index per row
    for i in range(len(Axis)): #prints row of data tied to index in header
        amount=int(graph[i,curr])
        percent=round((amount/mismatch)*100, len(str(mismatch))) #rounds % to sigfig based on number of mismatches
        output.write(f"{amount}({percent}%)\t") #data
    curr+=1
    output.write("\n") #move to next row

output.close

'''Graphing Graveyard''' #wanted to do a heatmap but ran out of time, will come back to this later, please ignore
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