#!/usr/bin/env python
#version 3.0
#Author: Farris Tedder

'''This module is a collection of bioinformatics functions written for the BGMP program at University of Oregon (Summer 2023).'''

###Constants###
DNA_bases = set('ATGCNatcgn')
RNA_bases = set('AUGCNaucgn')
compliment_DNA = {"A": "T", "T": "A", "C": "G", "G": "C", "U": "A", "N": "N"} #can be used with either RNA or DNA, produces DNA
compliment_RNA = {"A": "U", "T": "A", "C": "G", "G": "C", "U": "A", "N": "N"} #can be used with either RNA or DNA, produces RNA
###Functions###

def reverse_complement(line: str, RNAflag: bool=False) -> str:
    '''takes in a string of DNA or RNA, returns the reverse compliment'''
    newline = ""
    line=line.strip()
    line=line.upper()
    dict = (compliment_RNA if RNAflag else compliment_DNA)
    if validate_base_seq(line):
        for i in range(len(line)-1, -1, -1):
            newline+=dict[line[i]]
    return newline

def gc_content(DNA) -> float:
    '''Returns GC content of a DNA or RNA sequence as a decimal between 0 and 1.'''
    assert validate_base_seq(DNA), "String contains invalid characters - are you sure you used a DNA or RNA sequence?"
    DNA = DNA.upper()
    return (DNA.count("G")+DNA.count("C"))/len(DNA)

def validate_base_seq(seq: str, RNAflag: bool=False) -> bool: 
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag is set to true), Gs, Cs. False otherwise. Case insensitive.'''
    seq = seq.upper()
    return set(seq)<=(RNA_bases if RNAflag else DNA_bases)

def validate_DNA_seq(DNA) -> bool:
    '''Returns True if DNA is composed of A,C,G, and T only. Not case sensitive. False otherwise.'''
    DNA=DNA.upper()
    return (len(DNA) == DNA.count("A")+DNA.count("C")+DNA.count("G")+DNA.count("T")) #uses sum to see if any char in the string is not an A,T,C, or G

def convert_phred(phredchar) -> int:
    '''This funtion takes a valid askii charicter and returns the phred+33 score. Must be a string or int.'''
    if type(phredchar) != str:
        phredchar=str(phredchar) # sets everything to a string, for numbers
    return (ord(phredchar)-33)

def qual_score(phred_score: str) -> float:
    """Takes a string, return the average phred+33 quality score."""
    total = 0
    for letter in phred_score:
        score = convert_phred(letter)
        total = total+score
    return(total/len(phred_score))

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

def oneline_fasta(file: str, outputFileName: str):
    """Takes a FASTA file location and a path for output and returns to given output file a FASTA with all sequences as one line."""
    with open (file, "r") as FASTAfile, open (outputFileName, "w") as output: #parsing the raw data to put all sequences on one line
        totalLineNumber=0
        for line in FASTAfile:
            if line.startswith(">"): #if a header
                totalLineNumber+=1
                if totalLineNumber==1: #prints the first header without /n to avoid leading newline
                    output.write(line)
                else:
                    line = f"\n{line}"
                    output.write(line)
            else: #removes newlines from seq lines to make them 1 line
                output.write(line.strip("\n"))

def calc_median(sortedlist: list) -> float:
    """Takes in a sorted list and returns the median as a float."""
    mid=int(len(sortedlist)//2)
    if len(sortedlist)%2: #if length odd
        median = sortedlist[mid]
    else: #if length even
        median = ((sortedlist[mid]+sortedlist[mid-1])/2)
    return median

def init_list(length: int, value: float=0.0) -> list:
    '''This function takes an int and will populate a list of that length it with
    the value passed in "value". If no value is passed, initializes list with values of 0.0.'''
    i=0 #iterator
    finalList = []
    while i<length:
        finalList.append(value)
        i=i+1
    return finalList

###Functioning Tests###
if __name__ == "__main__":
    '''Base and DNA seq checkers'''
    assert validate_base_seq("AATAGAT") == True, "Validate base seq does not work on DNA"
    assert validate_base_seq("AAUAGAU", True) == True, "Validate base seq does not work on RNA"
    assert validate_base_seq("Hi there!") == False, "Validate base seq fails to recognize nonDNA"
    assert validate_base_seq("Hi there!", True) == False, "Validate base seq fails to recognize nonDNA"
    assert validate_DNA_seq("AATAGAT") == True, "DNA string not recognized"
    assert validate_DNA_seq("aatcga") == True, "DNA string not recognized"
    assert validate_DNA_seq("Hi there!") == False, "Non-DNA identified as DNA"
    assert reverse_complement("AGC") == "GCT", "incorrect complement"
    assert reverse_complement("ttca") == "TGAA", "incorrect complement"
    assert reverse_complement("AGC", True) == "GCU", "incorrect complement"
    '''quality score tester'''
    phred_score: str = "FFHHHHHJJJJIJIJJJIJJJJJJIIIJJJEHJJJJJJJIJIDGEHIJJFIGGGHFGHGFFF@EEDE@C??DDDDDDD@CDDDDBBDDDBDBDD@"
    assert qual_score("EEE") == 36
    assert qual_score("#I") == 21
    assert qual_score("EJ") == 38.5
    assert qual_score(phred_score) == 37.62105263157895, "wrong average phred score"
    '''GC Score'''
    assert gc_content("CGCGCGCCCG") == 1, "Wrong GC content score"
    assert gc_content("ATATATA") == 0, "Wrong GC content score"
    assert gc_content("AAAGGG") == 0.5, "Wrong GC content score"
    print("Passed DNA, RNA, and quality score varification tests")
    """Check that convert_phred returns the correct value for several different inputs"""
    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"
    '''quality score tester'''
    phred_score: str = "FFHHHHHJJJJIJIJJJIJJJJJJIIIJJJEHJJJJJJJIJIDGEHIJJFIGGGHFGHGFFF@EEDE@C??DDDDDDD@CDDDDBBDDDBDBDD@"
    assert qual_score("EEE") == 36
    assert qual_score("#I") == 21
    assert qual_score("EJ") == 38.5
    assert qual_score(phred_score) == 37.62105263157895, "wrong average phred score"
    print("Passed phred score conversion test")
    '''oneline_fasta tester'''
    import os
    inputpath="OneLineTestInputTester.txt"
    outputpath="OneLineTestOutputTester.txt"
    input=open(inputpath, "w+")
    line=""
    input.write(">This is the first header\n")
    input.write("ACTGACTGTACGT\n")
    input.write("ATCGATCGATCGA\n")
    input.write("> Header 2\n")
    input.write("ACTatgcatcATCGATC\n")
    input.write("TGCTAGCATC\n")
    input.write("TAGCAT\n")
    input.write("TAGCTAGA\n")
    input.close()
    oneline_fasta(inputpath, outputpath)
    output=open(outputpath, "r+")
    line=output.read()
    output.close()
    os.remove(outputpath)
    os.remove(inputpath)
    assert line == ">This is the first header\nACTGACTGTACGTATCGATCGATCGA\n> Header 2\nACTatgcatcATCGATCTGCTAGCATCTAGCATTAGCTAGA", "Output file incorrect"
    '''line checkers'''
    assert DNA_line(2) == True, "Is a DNA line"
    assert DNA_line(602) == True, "Is a DNA line"
    assert DNA_line(4) == False, "Not a DNA line"
    assert quality_score_line(4) == True, "Is a quality line"
    assert quality_score_line(240) == True, "Is a quality line"
    assert quality_score_line(2) == False, "Not a quality line"
    print("Passed line testers")
    '''Math'''
    testlist=[1,3,5,7,9]
    testlist2=[1,2,3,4,5,6,7,8,9,10]
    testlist3=[1,2,2,3]
    assert calc_median(testlist) == 5, "wrong median"
    assert calc_median(testlist2) == 5.5, "wrong median"
    assert calc_median(testlist3) == 2, "wrong median"    
    print("Passed math tests")
    '''List'''
    assert init_list(3) == [0.0, 0.0, 0.0], "incorrect list"
    assert init_list(3, 5) == [5.0, 5.0, 5.0], "incorrect list"
    print("List created correctly")
    '''reverse compliment'''
    assert reverse_compliment("ATG") == "CAT", "incorrect output"
    #broken, but not needed for current project. Fix later.
    '''one line from file
    from os import path
    oneline_fasta("./bioinfoonelinetestinput.txt", "./bioinfotestoutput.txt")
    assert path.exists("./bioinfotestoutput.txt") == True, "Output file not created"
    with open("./bioinfotestoutput.txt","r") as outputFile:
        line=str(outputFile.readline())
        assert line == "This is a multiline file that is composed of 5 lines that should be output as a single line with proper spacing. THISSHOULDALLBEONELINEWITHNOSPACESLIKEADNASEQ.", "Output file incorrect."
    '''
    print("Outputs created correctly")
