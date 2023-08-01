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
#### View file 1:
zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz | head
-Normal fastq file, 101 long seq, R1
-Header: @K00337:83:HJKJNBBXX:8:1101:1265:1191 1:N:0:1
#### View file 2:
zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz | head
only 8 char long seq, barcodes for R1
###### Record:
@K00337:83:HJKJNBBXX:8:1101:1265:1191 2:N:0:1
NCTTCGAC
+
\#AA\<FJJJ
#### View file 3:
zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz | head
only 8 char long seq, barcodes for R2
###### Record:
@K00337:83:HJKJNBBXX:8:1101:1265:1191 3:N:0:1
NTCGAAGA
+
\#AAAAJJF
#### View file 4:
zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz | head
-Normal fastq file, 101 long seq, R2
Header: @K00337:83:HJKJNBBXX:8:1101:1347:1191 4:N:0:1
###### Phred+33 encoded: # present in q-score, phred+64 has no #.
#### Count Ns
zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz | awk 'NR % 4 == 2' | tr -cd 'N' | wc -c
Run for files 2 and 3. Will not work if lower case n is present. All upper case in these file.
Index 1: 3976613
Index 2:  3329901
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
