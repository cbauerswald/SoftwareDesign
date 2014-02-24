# -*- coding: utf-8 -*-
"""
Created on Sun Feb  2 11:24:42 2014

@author: Cecelia Auerswald
"""

# you may find it useful to import these variables (although you are not required to use them)
from amino_acids import aa, codons
from random import shuffle

def collapse(L):
    """ Converts a list of strings to a string by concatenating all elements of the list """
    output = ""
    for s in L:
        output = output + s
    return output


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).
        
        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment
    """
    codon = ""
    proteins = ""
    for i in range(len(dna)): #looking at each of the letters
        codon += dna[i]
        if (i+1)%3==0: #if we have collected a full codon in the string codon
            for c in codons: 
                if codon in c: #check all possible codons for a match
                    proteins = proteins + aa[codons.index(c)] # This could be a +=, but this is fine too
            codon = "" #reset codon to check next set of three letters
    return proteins
            
        
    # Looks great - elegantly compact and well commented.

def coding_strand_to_AA_unit_tests():
    """ Unit tests for the coding_strand_to_AA function """

    print "input: TTT, expected output: F , actual output: "+coding_strand_to_AA("TTT")
    print "input: ACTGTCAAT, expected output: TVN , actual output: "+coding_strand_to_AA("ACTGTCAAT")
        
    # Nice unit tests throughout - dilligently done.

def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence
    
        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    """
    dna_list = list(dna)
    for i in range(len(dna)): #get the complement of each nucleotide 
        n=dna_list[i]
        if n == 'A':
            dna_list[i]='T'
        elif n == 'T':
            dna_list[i]='A'
        elif n=='C':
            dna_list[i] = 'G'
        elif n=='G':
            dna_list[i] = 'C'
    dna_list.reverse() #reverse the list of complementary nucleotides
    return ''.join(dna_list)

def get_reverse_complement_unit_tests():
    """ Unit tests for the get_complement function """
    print "input: ATGCCCGCTTT, expected output: AAAGCGGGCAT , actual output: "+ get_reverse_complement("ATGCCCGCTTT")
    print "input: CCGCGTTCA, expected output: TGAACGCGG , actual output: "+ get_reverse_complement("CCGCGTTCA")
    # YOUR IMPLEMENTATION HERE    

def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start codon and returns
        the sequence up to but not including the first in frame stop codon.  If there
        is no in frame stop codon, returns the whole string.
        
        dna: a DNA sequence
        returns: the open reading frame represented as a string
    """
    codon = ""
    for i in range(len(dna)): #looking at each of the letters
        codon += dna[i]
        if (i+1)%3==0: #if we have collected a full codon in the string codon
            if codon in ['TAG', 'TAA','TGA']: #check if codon in frame stop coding
                return dna[:i-2] #if so, return string up to that point
            codon=''
    return dna #if not stop codons found, return whole string
    # YOUR IMPLEMENTATION HERE

def rest_of_ORF_unit_tests():
    """ Unit tests for the rest_of_ORF function """
    print "input: ATGTGAA, expected output: ATG , actual output: "+ rest_of_ORF("ATGTGAA")
    print "input: ATGAGATAGG, expected output: ATGAGA, actual output: "+ rest_of_ORF("ATGAGATAGG")
    print "input: ATGGGGAGCTTG, expected output: ATGGGGAGCTTG, actual output: "+ rest_of_ORF("ATGGGGAGCTTG")
    
    
    # YOUR IMPLEMENTATION HERE
        
def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence and returns
        them as a list.  This function should only find ORFs that are in the default
        frame of the sequence (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.
        
        dna: a DNA sequence
        returns: a list of non-nested ORFs
    """
    codon=''
    allORFs=[]
    i=0
    while i<len(dna): #look through letters in dna sequence
        codon+=dna[i]
        if (i+1)%3==0:
            if codon=="ATG": #check for start codon
                allORFs.append(rest_of_ORF(dna[i-2:])) #if ORF found, add to list of ORFs
                i+=len(rest_of_ORF(dna[i-2:]))-3 #jump to end of that ORF in dna sequence
            codon = ''
        i+=1
    return allORFs
    # YOUR IMPLEMENTATION HERE     
     
def find_all_ORFs_oneframe_unit_tests():
    """ Unit tests for the find_all_ORFs_oneframe function """

    print "input: ATGCATGAATGTAGATAGATGTGCCC, expected output: ['ATGCATGAATGTAGA', 'ATGTGCCC'] , actual output: "+ str(find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC"))
    print "input: ATGTGCATGAGATAGGATGGGATGCTTG, expected output: ['ATGTGCATGAGA', 'ATGCTTG'], actual output: "+ str(find_all_ORFs_oneframe("ATGTGCATGAGATAGGATGGGATGCTTG"))
    
    # YOUR IMPLEMENTATION HERE

def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in all 3
        possible frames and returns them as a list.  By non-nested we mean that if an
        ORF occurs entirely within another ORF and they are both in the same frame,
        it should not be included in the returned list of ORFs.
        
        dna: a DNA sequence
        returns: a list of non-nested ORFs
    """
    return find_all_ORFs_oneframe(dna) + find_all_ORFs_oneframe(dna[1:])+find_all_ORFs_oneframe(dna[2:])
    
        
    # YOUR IMPLEMENTATION HERE

def find_all_ORFs_unit_tests():
    """ Unit tests for the find_all_ORFs function """
    
    # You are running the first test on a string other than your 
    # designated input string and as a result it is outputting values
    # that don't match your input values!

    print "input: ATGCATGAATGTAG, expected output: ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG'] , actual output: "+ str(find_all_ORFs("ATGCATGAATGTAGATAGATGTGCCC"))
    print "input: ATGTGCATGAGATAGGATGGGATGCTTG, expected output: ['ATGTGCATGAGA', 'ATGCTTG', 'ATGGGATGCTTG'], actual output: "+ str(find_all_ORFs("ATGTGCATGAGATAGGATGGGATGCTTG"))
    
    # YOUR IMPLEMENTATION HERE

def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.
        
        dna: a DNA sequence
        returns: a list of non-nested ORFs
    """
    
    return find_all_ORFs(dna) + find_all_ORFs(get_reverse_complement(dna))

    # Succintly done. Nice.

def find_all_ORFs_both_strands_unit_tests():
    """ Unit tests for the find_all_ORFs_both_strands function """
    
    # Same issue as in find_all_ORFs_unit_tests - make sure
    # that you are testing the same string you say you are!

    print "input: ATGCGAATGTAGCATCAAA, expected output: ['ATGCGAATG', 'ATGCTACATTCGCAT'] , actual output: "+ str(find_all_ORFs_both_strands("ATGCATGAATGTAGATAGATGTGCCC"))
    print "input: ATGTGCATGAGATAGGATGGGATGCTTG, expected output: ['ATGTGCATGAGA', 'ATGCTTG', 'ATGGGATGCTTG', 'ATGCACAT'], actual output: "+ str(find_all_ORFs_both_strands("ATGTGCATGAGATAGGATGGGATGCTTG"))
    
    # YOUR IMPLEMENTATION HERE

def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string"""
    ORFs = find_all_ORFs_both_strands(dna)
    longest_orf = ''
    for orf in ORFs:
        if len(orf)>len(longest_orf):
            longest_orf=orf
    return longest_orf

    # This works just fine, but for something like this it can 
    # be more compact to use the built-in max() function with len as a key.
    # I wouldn't typically mention something this nitpicky but for the
    # fact that your code is already very high quality.


def longest_ORF_unit_tests():
    """ Unit tests for the longest_ORF function """


    print "input: ATGCGAATGTAGCATCAAA, expected output: 'ATGCTACATTCGCAT', actual output: "+ str(longest_ORF("ATGCGAATGTAGCATCAAA"))
    
    # YOUR IMPLEMENTATION HERE

def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence
        
        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    longest_orf = ''
    dna_list = list(dna)
    for i in range(num_trials):
        shuffle(dna_list)
        orf = longest_ORF(collapse(dna_list))
        if len(orf)>len(longest_orf):
            longest_orf = orf
    return len(longest_orf)
        
    # Once again, this is a suggestion, not a critique: consider using list comprehensions
    # (look them up if you're not sure what they are, they're one of the coolest features of python).
   	# With that as a tool, you could go so far as to one-line this function with something like this:
    # return len(max([shuffle(list(dna)) for i in range(num_trials)]));

    # Also, this code is pretty readable, but it wouldn't hurt to comment here as diligently as at the beginning 

def gene_finder(dna, threshold):
    """ Returns the amino acid sequences coded by all genes that have an ORF
        larger than the specified threshold.
        
        dna: a DNA sequence
        threshold: the minimum length of the ORF for it to be considered a valid
                   gene.
        returns: a list of all amino acid sequences whose ORFs meet the minimum
                 length specified.
    """
    aa_seqs = []
    
    orfs = find_all_ORFs_both_strands(dna)
    orfs.sort(key=len, reverse=True)
    for orf in orfs:
        if len(orf)<=threshold:
            break
        aa_seqs.append(coding_strand_to_AA(orf))
    return aa_seqs
    
    
    
    
    
    
    
    # YOUR IMPLEMENTATION HERE