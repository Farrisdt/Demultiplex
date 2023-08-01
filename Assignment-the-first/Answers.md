# Assignment the First

## Part 1
1. Be sure to upload your Python script. Provide a link to it here:
[Data_Exploration.py](https://github.com/Farrisdt/Demultiplex/blob/a0f8dc9515f9e0c38dfec904f07e7214c4068ceb/Assignment-the-first/Data_Exploration.py)

| File name | label | Read length | Phred encoding |
|---|---|---|---|
| 1294_S1_L008_R1_001.fastq.gz | Read 1 | 101 | Phred+33 |
| 1294_S1_L008_R2_001.fastq.gz | Index 1 | 8 | Phred+33 |
| 1294_S1_L008_R3_001.fastq.gz | Index 2 | 8 | Phred+33 |
| 1294_S1_L008_R4_001.fastq.gz | Read 2 | 101 | Phred+33 |

2. Per-base NT distribution
    1. Use markdown to insert your 4 histograms here.
    ![Index 1](https://github.com/Farrisdt/Demultiplex/blob/3785180ffdac23f7f78845f2fbc7ad47a57da0f7/Assignment-the-first/Index1QualityScores.png)
    ![Index 2](https://github.com/Farrisdt/Demultiplex/blob/d9e14274947488bf69027a560544fe479d7c270d/Assignment-the-first/Index2QualityScores.png)
![Read 1](https://github.com/Farrisdt/Demultiplex/blob/d9e14274947488bf69027a560544fe479d7c270d/Assignment-the-first/Read1QualityScores.png)
![Read 2](https://github.com/Farrisdt/Demultiplex/blob/d9e14274947488bf69027a560544fe479d7c270d/Assignment-the-first/Read2QualityScores.png)
    ii. What is a good quality score cutoff for index reads and biological read pairs to utilize for sample identification and downstream analysis, respectively? Justify your answer.
    
    30 is a good quality score cut off because it is just below the lowest averages per base. Most of the average quality scores are around the same value so a higher quality score is doable. Using an average qscore for the entire read will allow us to not filter out reads that only have one or two low qscore bases. Since we are already filtering out all sequences with N, unknown bases will not bring down the average read score. Checking per base is an option, but since the likelihood of a false positive from a single base substitution is low, and since the qscore would lower below the threshold if half of the bases are of low quality, I think taking the average of the entire set is appropriate. Something more complex, such as the average of the lowest 4 scores per read, would give even more refinement but it does not feel necessary for this project. I was thinking 30 would be the minimum so a 29 would be filtered out. If the data is looking off raising this value to a 33 or even 35 may be appropriate, but I think it is best to start conservatively and remove data as needed. 
       
    iii. How many indexes have undetermined (N) base calls? (Utilize your command line tool knowledge. Submit the command(s) you used. CHALLENGE: use a one-line command)

   Index 1: 3976613

   Index 2: 3329901
   
```
/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz | awk 'NR % 4 == 2' | tr -cd 'N' | wc -c
```

## Part 2
1. Define the problem
2. Describe output
3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [>=6 expected output FASTQ files](../TEST-output_FASTQ).
4. Pseudocode
5. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement
