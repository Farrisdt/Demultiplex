## Class Notes
Leslie: De-multiplexing:

Text file with indexes on talapas. 4 fastq files in bgmp/shared/2017_sequencing. DO NOT COPY, DO NOT UNZIP, ONLY LOOP THROUGH ONCE, look into gzip for python.

Files 1+4 are reads, 2+3 are barcodes. Reads have 8nt twin indexs that are ‘duel matched’. reads 3 and 4 must be reverse complemented due to the bridging amplification (to match the index list, 5’first). Unknown nucleotides (N) can still be found out if they match enough of a single known index, reversed barcode can also be used to suss out missing nt using it’s match and it’s matches quality score.

For the sequence, since it is at opposite sides of the strand they will not be reverse complimented.

If confused, look up “illumina sequencing by synthesis”.

Key words: illumina, sequencing by synthesis, bridge amplification, index/barcode, duel matched indexes, de-multiplexing

               Part 1: Quality score distribution per nucleotide:

Plot base position (x) vs average quality score (y). gather read length, phred encoding, and label. Use running sum strategy

Submit lab notebook with all bash scripts (yes even wc -l)

Header: @K00337:83:HJKJNBBXX:8:1101:1265:1191 1:N:0:1

               Green: read number (1-4)

Sudo code discussion:

Need a R1 and R2 file for each properly matched index (minimum 2 for if only 1 index was matched), as well as 2 for hopped and unknown (6 total files per index)

Assignment 2: make issues on 3 classmates code on github
# Assignment One
## Part 1:
Run in bgmp_py311
installed matplotlib to environment, most recent version as of 7/31/2023

	conda install matplotlib
### Bash:
	#!/bin/bash
	
	#SBATCH --account=bgmp
	
	#SBATCH --partition=compute
	
	conda activate bgmp_py311
	
	  
	
	./Data_Exploration.py -f "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz" -o "Index1QualityScores" -t 2 -l 8
	
	./Data_Exploration.py -f "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz" -o "Index2QualityScores" -t 3 -l 8
	
	./Data_Exploration.py -f "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz" -o "Read1QualityScores" -t 1 -l 101
	
	./Data_Exploration.py -f "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz" -o "Read2QualityScores" -t 4 -l 101
### Data Exploration:
#### To get char number per sequence: 
zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz | head -n 4 | tail -n 1 | wc -m
-102
-(101 non-new line chars, 101 nucleotides)
8 char for barcodes, given information confirmed by hand counting data below.
#### To get lines in file:
zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz | wc -l
1452986940
###### Double check with barcode:
Using the r1 barcodes checks both that barcodes have the same number of lines as reads and that read 1 and 2 are the same.
zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz | wc -l
1452986940
#### View Files
##### View file 1:
zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz | head
-Normal fastq file, 101 long seq, R1
-Header: @K00337:83:HJKJNBBXX:8:1101:1265:1191 1:N:0:1
##### View file 2:
zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz | head
only 8 char long seq, barcodes for R1
###### Record:
@K00337:83:HJKJNBBXX:8:1101:1265:1191 2:N:0:1
NCTTCGAC
+
\#AA\<FJJJ
##### View file 3:
zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz | head
only 8 char long seq, barcodes for R2
###### Record:
@K00337:83:HJKJNBBXX:8:1101:1265:1191 3:N:0:1
NTCGAAGA
+
\#AAAAJJF
##### View file 4:
zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz | head
-Normal fastq file, 101 long seq, R2
Header: @K00337:83:HJKJNBBXX:8:1101:1347:1191 4:N:0:1
###### Phred+33 encoded: # present in q-score, phred+64 has no #.
#### Count Ns
zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz | awk 'NR % 4 == 2' | grep-c 'N'
Run for files 2 and 3. Will not work if lower case n is present. All upper case in these file.
Index 1: 3976613
Index 2:  3328051
#### Table

| File name |  label  | Read length | Phred encoding |
|---|---|---|---|
| 1294_S1_L008_R1_001.fastq.gz | Read 1 | 101 | Phred+33 |
| 1294_S1_L008_R2_001.fastq.gz | Index 1 | 8 | Phred+33 |
| 1294_S1_L008_R3_001.fastq.gz | Index 2 | 8 | Phred+33 |
| 1294_S1_L008_R4_001.fastq.gz | Read 2 | 101 | Phred+33 |

## Part 2:
#### Testers:
 Used http://www.faculty.ucr.edu/~mmaduro/random.htm to generate random DNA sequences and http://www.unit-conversion.info/texttools/random-string-generator/ to generate random quality scores.

# Assignment Three
Assignment 2 was to give feedback to others, no notes needed.
## Discussion Notes:
% of reads and raw count from each sample
amount of index hopping
amount of index hopping between known indexes, create a table/heatmap for all cross sections (like a multiplication table) be sure to include both index1vsindex4 and index4vsindex1
-talk about no way of knowing if both strands are index hopped to the same wrong index
## Bash
### Code
	#!/bin/bash
	/#SBATCH --account=bgmp
	/#SBATCH --partition=compute
	
	conda activate bgmp_py311 /#activate environment

	/usr/bin/time -v ./code_code.py -r1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz -r2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz -i1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz -i2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz

	/usr/bin/time -v ./code_code.py -q 35 -r1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz -r2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz -i1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz -i2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz
	
	/usr/bin/time -v ./code_code.py -q 25 -r1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz -r2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz -i1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz -i2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz
	
	/usr/bin/time -v ./code_code.py -q 20 -r1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz -r2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz -i1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz -i2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz
	
	/usr/bin/time -v ./code_code.py -q 0 -r1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz -r2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz -i1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz -i2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz

Submitted batch job 26469
	Slurm output location:
	
	/projects/bgmp/ftedder/bioinfo/Bi622/Demultiplex/Assignment-the-third/slurm-26469.out 
### Notes
Ran quality scores 0, 20, 25, 30, and 35.
The desired quality score can change depending on what the data is being used for. Since I am currently unsure of what is needed down stream I offered a selection of results. 25-30 is generally the best cut off to ensure primarily high quality data. 35 is highly stringent but could be needed for certain projects where error is especially troublesome. 20 is low, but has more data if error is not a worry. 0 is added mostly in case anyone wants to see everything or as a good comparison for uncertainty due to 'Ns' vs low quality scores.
# General Notes:
Installed/updated matplotlib to env bgmp_py311:

	The following packages will be downloaded:

    package                    |            build
    ---------------------------|-----------------
    openssl-3.1.2              |       hd590300_0         2.5 MB  conda-forge
    ------------------------------------------------------------
                                           Total:         2.5 MB

	The following packages will be UPDATED:

	  openssl                                  3.1.1-hd590300_1 --> 3.1.2-hd590300_0 

Updated bioinfo.py to version 3. Added compliment dictionaries for RNA and DNA, as well as a reverse compliment function that can output either RNA or DNA.

## Python
For Python: code_code.py (I don't want to change the name, I've bonded with it.)
args: (need to know if running)
	
	("-r1", "--read1", help="read 1 file name", type=str, required=True)
    ("-r2", "--read2", help="read 2 file name", type=str, required=True)
    ("-i1", "--index1", help="index 1 file name", type=str, required=True)
	("-i2", "--index2", help="index 2 file name", type=str, required=True)        
    ("-i", "--indexL", help="index length", default=8, type=int, required=False)
    ("-q", "--qualcut", help="quality score cutoff", default=30, type=int, required=False)
    ("-l", "--lineNumber", help="number of lines in the file", default=1452986940, type=int, required=False)

There is a graphics graveyard at the bottom of the code for when the heatmap was going to be colored instead of just a chart. If this ends up being code I use I will fix and implement it, but as of now it's not needed. I believe the chart should be easily heat-mappable in the future with the way I've set it up. 
###### Heatmap idea:
Just make an array of color values while making the chart. So as looping have a for statement that just grabs the value (percent or raw or both) and then if <10 make blue, 10-20 make green, etc.
