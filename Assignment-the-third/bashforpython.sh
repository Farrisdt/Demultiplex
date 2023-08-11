#!/bin/bash
#SBATCH --account=bgmp
#SBATCH --partition=compute

conda activate bgmp_py311

/usr/bin/time -v ./code_code.py -r1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz -r2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz -i1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz -i2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz
/usr/bin/time -v ./code_code.py -q 35 -r1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz -r2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz -i1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz -i2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz
/usr/bin/time -v ./code_code.py -q 25 -r1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz -r2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz -i1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz -i2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz
/usr/bin/time -v ./code_code.py -q 20 -r1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz -r2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz -i1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz -i2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz
/usr/bin/time -v ./code_code.py -q 0 -r1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz -r2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz -i1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz -i2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz

#For testing
#./code_code.py -l 20 -r1 ../TEST-input_FASTQ/read1test.fastq -r2 ../TEST-input_FASTQ/read2test.fastq -i1 ../TEST-input_FASTQ/index1test.fastq -i2 ../TEST-input_FASTQ/index2test.fastq