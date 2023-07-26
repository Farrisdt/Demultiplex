#!/usr/bin/env python
#version 2.0: most current verion 7/26/2023
#Author: Farris Tedder, for Bi621

'''This module is a collection of bioinformatics functions written for Bi621 at University of Oregon (Summer 2023).'''

###Constants###
DNA_bases = set('ATGCNatcgn')
RNA_bases = set('AUGCNaucgn')

###Functions###
def gc_content(DNA):
    '''Returns GC content of a DNA or RNA sequence as a decimal between 0 and 1.'''
    assert validate_base_seq(DNA), "String contains invalid characters - are you sure you used a DNA or RNA sequence?"
    
    DNA = DNA.upper()
    return (DNA.count("G")+DNA.count("C"))/len(DNA)

def validate_base_seq(seq: str, RNAflag: bool=False) -> bool: 
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    seq = seq.upper()
    return set(seq)<=(RNA_bases if RNAflag else DNA_bases)

def validate_DNA_seq(DNA) -> bool:
    '''Returns True if DNA is composed of A,C,G, and T only. Not case sensitive. False otherwise.'''
    '''DNA =DNA.upper()
    length = len(DNA)
    array = DNA.count("A")+DNA.count("C")+DNA.count("G")+DNA.count("T")'''
    DNA=DNA.upper()
    # if (DNA == [ATGC]*){return True}
    return (len(DNA) == DNA.count("A")+DNA.count("C")+DNA.count("G")+DNA.count("T"))
    # use 'pass' to bipass incommplete fxn errors while coding, will return nothing (none)

def convert_phred(phredchar) -> int:
    '''This funtion takes a valid askii charicter and returns the phred33 score. Must be a string or int.'''
    if type(phredchar) != str:
        phredchar=str(phredchar)
    return (ord(phredchar)-33)

def qual_score(phred_score: str) -> float:
    """Takes a string, return the average phred33 quality score."""
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

def calc_median (sortedlist: list) -> float:
    """Takes in a sorted list and returns a median float."""
    mid=int(len(sortedlist)//2)
    if len(sortedlist)%2: #if length odd
        median = sortedlist[mid]
    else: #if length even
        median = ((sortedlist[mid]+sortedlist[mid-1])/2)
    return median

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
    from os import path
    oneline_fasta("./bioinfoonelinetestinput.txt", "./bioinfotestoutput.txt")
    assert path.exists("./bioinfotestoutput.txt") == True, "Output file not created"
    with open("./bioinfotestoutput.txt","r") as outputFile:
        line=str(outputFile.readline())
        assert line == "This is a multiline file that is composed of 5 lines that should be output as a single line with proper spacing. THISSHOULDALLBEONELINEWITHNOSPACESLIKEADNASEQ.", "Output file incorrect."
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