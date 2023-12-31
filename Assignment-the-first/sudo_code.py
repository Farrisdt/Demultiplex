#!/usr/bin/env python

Problem: Need to know which records are index swapped/have undetermined indexes. 
Output: 52 files. One for each read for all pairs with an index not found in the list or of too low quality. One for each read for pairs who indexes do not match. One each for read 1 and read 2 per index. (24 indexs X 2 + 2 unknown + 2 mismatch = 52). Also the number of read-pairs in matched, unmatched, and unknown files.

def demultiplex(r1: str, r2: str, r3: str, r4: str):
    '''Takes in the paths, as strings, to the 4 file outputs of illumina seqeuncing. Outputs 52 files, sorted by read number and index including a mismatch file and an unknown file.'''
    create a dictionary of the given indexs. Can be used to hold the values of the number of read pairs written to each file. Each index could have a 3 spot list, [match, mismatch, unknown]
    open all 52 output files, not with 'with open', use a loop instead. Can loop through index dict keys and use them are variable names, then open read 1 and read 2 files for each index. Also open 1 unmatched and 1 unknown file outside the loop.
    #constants/variables
    i= 1452986940 #numbner of lines in file, found in data exploration
    quality cutoff score
    2 0-set int variabled used for checking quality score (averageQualityScore1 and averageQualityScore2)

    open all 4 data files
        loop through each line:
            read in line from each file, set to variable (such as read1, index1, etc.)
            if a header line: #could find @ line starter. may also make bioinfo script to check which data line we are on based on the line number, return true if header, like the ones I already have
                header1 = read1line
                header2 = read2line
            if a sequence line: #already have coded from Bi621
                reverse index2 and find the complement. Make a bioinfo script to do this automatically. Turn A>T etc and N>N, return new compliment string. make input upper case and error if not DNA (any char but ACTGN).
                append "index1-index2" to both headers
            if quality score line: #already in bioinfo
                reset qualityscoreaverage to 0 for index 1 and 2.
                sum the quality score of all char in the sequence into the variable for 1 and 2
            if both index1 and 2 are found in indexs dictionary:
                if both indexs quality score are above quality cutoff:
                    if indexs match:
                        write to read1 and read2 matched files
                        add 1 to the first index of the list in the index dict.
                    else:
                        write to mismatch file
                        add 1 to the second index of the list in the index dict
                else:
                    write to unknown file (due to failing quality score)
                    add 1 to the third index of the list in the index dict
            else:
                write to unknown file (due to unknown index)
                add 1 to the third index of the list in the index dict

    report out the variables stored in the index dict. (the number of read-pairs for each type of file)
    close all 52 files by using an edited version of the same loop from the begining.

def reverse_compliment_barcode(seq: str) -> str:
    '''Takes a index 2 from an illumina sequence and creates it's reverse compliment. Returns a string that should match index 1 exactly.'''
    return reversed reverse compliment
Input: AACGTNA
Expected output: TNACGTT

see tests at bottom for expected outputs for line checkers.
def header_line(i) -> bool:
    """checks to see if line from a FASTQ file is a header, returns bool. Takes the line number, 1 based index."""
    if i%4 == 1:
        return True
    else:
        return False

From bioinfo.py:
def DNA_line(i) -> bool:
    """checks to see if line from a FASTQ file is a DNA seqence, returns bool. Takes the line number, 1 based index."""
    if i%4 == 2:
        return True
    else:
        return False

def quality_score_line(i) -> bool:
    """checks to see if line from a FASTQ file is a phred quality score, returns bool. Takes the line number, 1 based index."""
    if i%4 == 0:
        return True
    else:
        return False
From bioinfo.py testers:

'''line checkers'''
    assert DNA_line(2) == True, "Is a DNA line"
    assert DNA_line(602) == True, "Is a DNA line"
    assert DNA_line(4) == False, "Not a DNA line"
    assert quality_score_line(4) == True, "Is a quality line"
    assert quality_score_line(240) == True, "Is a quality line"
    assert quality_score_line(2) == False, "Not a quality line"
    print("Passed line testers")

New possible tester:
    assert header_line(1) == True, "Is a header line"
    assert header_line(261) == True, "Is a header line"
    assert header_line(4) == False, "Not a header line"
