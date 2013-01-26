# encoding=utf8
import sys
import gzip
import os.path as path
from itertools import islice

"""
play with Rosalind http://rosalind.info
"""

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def count_DNA(DNA_string):
    """
Counting DNA Nucleotides
http://rosalind.info/problems/dna/

Problem

A string is simply an ordered collection of symbols selected from some alphabet and formed into a word; the length of a string is the number of symbols that it contains.

An example of a length 21 DNA string (whose alphabet contains the symbols 'A', 'C', 'G', and 'T') is "ATGCTTCAGAAAGGTCTTACG."

Given: A DNA string s of length at most 1000 nt.

Return: Four integers (separated by spaces) counting the respective number of times that the symbols 'A', 'C', 'G', and 'T' occur in s.
Sample Dataset

AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC

Sample Output

20 12 17 21

>>> count_DNA("AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC")
(20, 12, 17, 21)
"""    

    return DNA_string.count("A"), DNA_string.count("C"), DNA_string.count("G"), DNA_string.count("T")

#-------------------------------------------------------------------------------
def transcribe(DNA_string):
    """
Transcribing DNA into RNA
http://rosalind.info/problems/rna/

Problem

An RNA string is a string formed from the alphabet containing 'A', 'C', 'G', and 'U'.

Given a DNA string t corresponding to a coding strand, its transcribed RNA string u is formed by replacing all occurrences of 'T' in t with 'U' in u.

Given: A DNA string t having length at most 1000 nt.

Return: The transcribed RNA string of t.
Sample Dataset

GATGGAACTTGACTACGTAAATT

Sample Output

GAUGGAACUUGACUACGUAAAUU

>>> transcribe("GATGGAACTTGACTACGTAAATT")
'GAUGGAACUUGACUACGUAAAUU'
"""    

    return DNA_string.replace("T", "U")

#-------------------------------------------------------------------------------
import string
complement_trans = string.maketrans("ATGC", "TACG")

def reverse_complement_DNA(DNA_string):
    """
Complementing a Strand of DNA
http://rosalind.info/problems/revc/

Problem

In DNA strings, symbols 'A' and 'T' are complements of each other, as are 'C' and 'G'.

The reverse complement of a DNA string s is the string sc formed by reversing the symbols of s, then taking the complement of each symbol (e.g., the reverse complement of "GTCA" is "TGAC").

Given: A DNA string s of length at most 1000 bp.

Return: The reverse complement sc of s.
Sample Dataset

AAAACCCGGT

Sample Output

ACCGGGTTTT

>>> reverse_complement_DNA("AAAACCCGGT")
'ACCGGGTTTT'
"""

    return DNA_string.translate(complement_trans)[::-1]

#-------------------------------------------------------------------------------
def max_gc_content(fasta):
    r"""
Computing GC Content
http://rosalind.info/problems/gc/

    Identifying Unknown DNA Quicklyclick to collapse

    A quick method used by early computer software to determine the language of a given piece of text was to analyze the frequency with which each letter appeared in the text. This strategy was used because each language tends to exhibit its own letter frequencies, and as long as the text under consideration is long enough, software will correctly recognize the language quickly and with a very low error rate. See Figure 1 for a table compiling English letter frequencies.

    You may ask: what in the world does this linguistic problem have to do with biology? Although two members of the same species will have different genomes, they still share the vast percentage of their DNA; notably, 99.9% of the 3.2 billion base pairs in a human genome are common to almost all humans (i.e., excluding people having major genetic defects). For this reason, biologists will speak of the human genome, meaning an average-case genome derived from a collection of individuals. Such an average case genome can be assembled for any species, a challenge that we will soon discuss.

    The biological analog of identifying unknown text arises when researchers encounter a molecule of DNA deriving from an unknown species. Because of the base pairing relations of the two DNA strands, cytosine and guanine will always appear in equal amounts in a double-stranded DNA molecule. Thus, to analyze the symbol frequencies of DNA for comparison against a database, we compute the molecule's GC-content, or the percentage of its bases that are either cytosine or guanine.

    In practice, the GC-content of most eukaryotic genomes hovers around 50%. However, because genomes are so long, we may be able to distinguish species based off very small discrepancies in GC-content; furthermore, most prokaryotes have a GC-content significantly higher than 50%, so that GC-content can be used to quickly differentiate many prokaryotes and eukaryotes by using relatively small DNA samples.

Problem

The GC-content of a DNA string is given by the percentage of symbols in the string that are 'C' or 'G'. For example, the GC-content of "AGCTATAG" is 37.5%. Note that the reverse complement of any DNA string has the same GC-content.

DNA strings must be labeled when they are consolidated into a database. A commonly used method of string labeling is called FASTA format. In this format, the string is introduced by a line that begins with '>', followed by some labeling information. Subsequent lines contain the string itself; the first line to begin with '>' indicates the label of the next string.

In Rosalind's implementation, a string in FASTA format will be labeled by the ID "Rosalind_xxxx", where "xxxx" denotes a four-digit code between 0000 and 9999.

Given: At most 10 DNA strings in FASTA format (of length at most 1 kbp each).

Return: The ID of the string having the highest GC-content, followed by the GC-content of that string. Rosalind allows for a default error of 0.001 in all decimal answers unless otherwise stated; please see the note on absolute error below.
Sample Dataset

>Rosalind_6404
CCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCC
TCCCACTAATAATTCTGAGG
>Rosalind_5959
CCATCGGTAGCGCATCCTTAGTCCAATTAAGTCCCTATCCAGGCGCTCCGCCGAAGGTCT
ATATCCATTTGTCAGCAGACACGC
>Rosalind_0808
CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGAC
TGGGAACCTGCGGGCAGTAGGTGGAAT

Sample Output

Rosalind_0808
60.919540%

    Note on Absolute Errorclick to collapse

    We say that a number x is within an absolute error of y to a correct solution if x is within y of the correct solution. For example, if an exact solution is 6.157892, then for x to be within an absolute error of 0.001, we must have that |x−6.157892|<0.001, or 6.156892<x<6.158892.

    Error bounding is a vital practical tool because of the inherent round-off error in representing decimals in a computer, where only a finite number of decimal places are allotted to any number. After being compounded over a number of operations, this round-off error can become evident. As a result, rather than testing whether two numbers are equal with x=z, you may wish to simply verify that |x−z| is very small.

    The mathematical field of numerical analysis is devoted to rigorously studying the nature of computational approximation.

    >>> s = '>Rosalind_6404\nCCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCC\nTCCCACTAATAATTCTGAGG\n>Rosalind_5959\nCCATCGGTAGCGCATCCTTAGTCCAATTAAGTCCCTATCCAGGCGCTCCGCCGAAGGTCT\nATATCCATTTGTCAGCAGACACGC\n>Rosalind_0808\nCCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGAC\nTGGGAACCTGCGGGCAGTAGGTGGAAT\n'
    >>> max_gc_content(s)
    Rosalind_0808
    60.919540%
"""    

    gc = [(seq_name, (seq.count("G") + seq.count("C")) * 100.0 / len(seq)) for seq_name, seq in fasta_iter(fasta)]
    seq_name, value = max(gc, key=lambda ln: ln[1])
    result = "%s\n%.6f%%\n" % (seq_name, value)
    sys.stdout.write(result)
    sys.stdout.flush()

def fasta_iter(fa, buffsize=100000):
    r"""
    iter over a fasta file or file-like object or string.

    input: fa could be a file object, a filename or a string of fasta records
    return: record of fasta: ("sequence_name", "sequence")

    >>> s = '>Rosalind_6404\nCCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCC\nTCCCACTAATAATTCTGAGG\n>Rosalind_5959\nCCATCGGTAGCGCATCCTTAGTCCAATTAAGTCCCTATCCAGGCGCTCCGCCGAAGGTCT\nATATCCATTTGTCAGCAGACACGC\n>Rosalind_0808\nCCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGAC\nTGGGAACCTGCGGGCAGTAGGTGGAAT\n'
    >>> list(fasta_iter(s))
    [('Rosalind_6404', 'CCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCCTCCCACTAATAATTCTGAGG'), ('Rosalind_5959', 'CCATCGGTAGCGCATCCTTAGTCCAATTAAGTCCCTATCCAGGCGCTCCGCCGAAGGTCTATATCCATTTGTCAGCAGACACGC'), ('Rosalind_0808', 'CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGACTGGGAACCTGCGGGCAGTAGGTGGAAT')]
    """

    basestring_type = basestring if sys.version_info[0] == 2 else str

    if isinstance(fa, basestring_type):  # if fa is filename or fasta string 
        if path.exists(fa):              # fa is a file on disk
            if fa.endswith(".gz"): fobj = gzip.open(fa, 'rb')
            else: fobj = open(fa)
        elif fa == '-':                  # fa is standard input
            fobj = sys.stdin
            buffsize = 10000   # send output earlier if the input is sys.stdin
        elif fa[0] == ">" and "\n" in fa:# fa is a fasta string
            fobj = iter(fa.splitlines())
        else:
            raise Exception("Don't recognize the format. Valid formats includes: string, list/tuple of string, file-like object.")
    # if (fa is list or tuple of string), or (fa is file-like object)
    # in both cases, we can use "for line in fa" to iter over fa
    elif (hasattr(fa, '__getitem__') and callable(fa.__getitem__)) or (hasattr(fa, 'readline') and callable(fa.readline)):
        fobj = fa
    else:
        raise Exception("Don't recognize the format. Valid formats includes: string, list/tuple of string, file-like object.")

    chunk = []
    while True:
        new_data = list(islice(fobj, buffsize))
        if not new_data: break
        chunk.extend(new_data)
        idx = [i for i, ln in enumerate(chunk) if ln[0] == '>']
        for i, j in enumerate(idx[:-1]):
            yield (chunk[j][1:].rstrip('\n'), "".join(chunk[j+1:idx[i+1]]).replace("\n", ""))
        chunk = chunk[idx[-1]:]

    if chunk:
        yield (chunk[0][1:].rstrip('\n'), "".join(chunk[1:]).replace("\n", ""))

#-------------------------------------------------------------------------------
def count_point_mutation(s, t=None):
    """
Counting Point Mutations
http://rosalind.info/problems/hamm/

Problem

Given two strings s and t of equal length, the Hamming distance between s and t, denoted dH(s,t), is the number of corresponding symbols that differ in s and t. See Figure 2.

Given: Two DNA strings s and t of equal length (not exceeding 1 kbp).

Return: The Hamming distance dH(s,t).
Sample Dataset

GAGCCTACTAACGGGAT
CATCGTAATGACGGCCT

Sample Output

7

    >>> count_point_mutation('GAGCCTACTAACGGGAT', 'CATCGTAATGACGGCCT')
    7
"""    

    if t is None:
        s, t = s.split()

    if len(s) != len(t):
        raise Exception("two strings should be of the same length!")

    return sum(1 for i in range(len(s)) if s[i] != t[i])

#-------------------------------------------------------------------------------
def find_motif(seq, motif):
    """
Finding a Motif in DNA
http://rosalind.info/problems/subs/

    Combing Through the Haystackclick to collapse

    Finding the same interval of DNA in the genomes of two different organisms (often taken from different species) is highly suggestive that the interval has the same function in both organisms.

    We define a motif as such a commonly shared interval of DNA. A common task in molecular biology is to search an organism's genome for a known motif.

    The situation is complicated by the fact that genomes are riddled with intervals of DNA that occur multiple times (possibly with slight modifications), called repeats. These repeats occur far more often than would be dictated by random chance, indicating that genomes are anything but random and in fact illustrate that the language of DNA must be very powerful (compare with the frequent reuse of common words in any human language).

    The most common repeat in humans is the Alu repeat, which is approximately 300 bp long and recurs around a million times throughout every human genome. However, Alu has not been found to serve a positive purpose, and appears in fact to be parasitic: when a new Alu repeat is inserted into a genome, it frequently causes genetic disorders.

Problem

Given two strings s and t, t is a substring of s if t is contained as a contiguous collection of symbols in s (as a result, t must be no longer than s).

The position of a symbol in a string is the total number of symbols found to its left, including itself (e.g., the positions of all occurrences of 'U' in "AUGCUUCAGAAAGGUCUUACG" are 2, 5, 6, 15, 17, and 18). The symbol at position i of s is denoted by s[i].

A substring of s can be represented as s[j:k], where j and k represent the starting and ending positions of the substring in s; for example, if s = "AUGCUUCAGAAAGGUCUUACG", then s[2:5] = "UGCU".

The location of a substring s[j:k] is its beginning position j; note that t will have multiple locations in s if it occurs more than once as a substring of s (see the Sample below).

Given: Two DNA strings s and t (each of length at most 1 kbp).

Return: All locations of t as a substring of s.
Sample Dataset

GATATATGCATATACTT
ATAT

Sample Output

2 4 10

    >>> find_motif('GATATATGCATATACTT', 'ATAT')
    [2, 4, 10]
"""

    import re
    
    return [match.start()+1 for match in re.finditer(r'(?=%s)'%re.escape(motif), seq)]

#-------------------------------------------------------------------------------
def cons(s):
    """
Consensus and Profile
http://rosalind.info/problems/cons/

    Finding a Most Likely Common Ancestorclick to collapse

    In “Counting Point Mutations”, we calculated the minimum number of symbol mismatches between two strings of equal length to model the problem of finding the minimum number of point mutations occurring on the evolutionary path between two homologous strands of DNA. If we instead have several homologous strands that we wish to analyze simultaneously, then the natural problem is to find an average-case strand to represent the most likely common ancestor of the given strands.

Problem

A matrix is a rectangular table of values divided into rows and columns. An m×n matrix has m rows and n columns. Given a matrix A, we write Ai,j to indicate the value found at the intersection of row i and column j.

Say that we have a collection of DNA strings, all having the same length n. Their profile matrix is a 4×n matrix P in which P1,j represents the number of times that 'A' occurs in the jth position of one of the strings, P2,j represents the number of times that C occurs in the jth position, and so on (see below).

A consensus string c is a string of length n formed from our collection by taking the most common symbol at each position; the jth symbol of c therefore corresponds to the symbol having the maximum value in the j-th column of the profile matrix. Of course, there may be more than one most common symbol, leading to multiple possible consensus strings.

           	    A T C C A G C T
           	    G G G C A A C T
           	    A T G G A T C T
DNA Strings	    A A G C A A C C
           	    T T G G A A C T
           	    A T G C C A T T
           	    A T G G C A C T

           	A   5 1 0 0 5 5 0 0
    Profile	C   0 0 1 4 2 0 6 1
           	G   1 1 6 3 0 1 0 0
           	T   1 5 0 0 0 1 1 6
  Consensus	    A T G C A A C T

Given: A collection of at most 10 DNA strings of equal length (at most 1 kbp).

Return: A consensus string and profile matrix for the collection. (If several possible consensus strings exist, then you may return any one of them.)
Sample Dataset

ATCCAGCT
GGGCAACT
ATGGATCT
AAGCAACC
TTGGAACT
ATGCCATT
ATGGCACT

Sample Output

ATGCAACT
A: 5 1 0 0 5 5 0 0
C: 0 0 1 4 2 0 6 1
G: 1 1 6 3 0 1 0 0
T: 1 5 0 0 0 1 1 6

    >>> s = '''
    ... ATCCAGCT
    ... GGGCAACT
    ... ATGGATCT
    ... AAGCAACC
    ... TTGGAACT
    ... ATGCCATT
    ... ATGGCACT'''
    >>> cons(s)
    ATGCAACT
    A: 5 1 0 0 5 5 0 0
    C: 0 0 1 4 2 0 6 1
    G: 1 1 6 3 0 1 0 0
    T: 1 5 0 0 0 1 1 6
"""

    if isinstance(s, basestring):
        d = s.split()
    elif isinstance(s, (list, tuple)):
        d = s
    else:
        raise Exception("invalid input")

    ACGT = 'ACGT'

    from collections import defaultdict
    mat = defaultdict(list)
    consensus = []

    mat = zip(*[[col.count(C) for C in ACGT] for col in zip(*d)])
    consensus = "".join([ACGT[col.index(max(col))] for col in zip(*mat)])

    print consensus
    for C, ln in zip(ACGT, mat):
        print "%s: %s" % (C, " ".join(map(str, ln)))

    f = open('cons.txt', 'w')
    f.write("%s\n" % consensus)
    for C, ln in zip(ACGT, mat):
        f.write("%s: %s\n" % (C, " ".join(map(str, ln))))
    f.close()

#-------------------------------------------------------------------------------
RNA2protein = {
    "UUU": "F",      "CUU": "L",      "AUU": "I",      "GUU": "V",
    "UUC": "F",      "CUC": "L",      "AUC": "I",      "GUC": "V",
    "UUA": "L",      "CUA": "L",      "AUA": "I",      "GUA": "V",
    "UUG": "L",      "CUG": "L",      "AUG": "M",      "GUG": "V",
    "UCU": "S",      "CCU": "P",      "ACU": "T",      "GCU": "A",
    "UCC": "S",      "CCC": "P",      "ACC": "T",      "GCC": "A",
    "UCA": "S",      "CCA": "P",      "ACA": "T",      "GCA": "A",
    "UCG": "S",      "CCG": "P",      "ACG": "T",      "GCG": "A",
    "UAU": "Y",      "CAU": "H",      "AAU": "N",      "GAU": "D",
    "UAC": "Y",      "CAC": "H",      "AAC": "N",      "GAC": "D",
    "UAA": "Stop",   "CAA": "Q",      "AAA": "K",      "GAA": "E",
    "UAG": "Stop",   "CAG": "Q",      "AAG": "K",      "GAG": "E",
    "UGU": "C",      "CGU": "R",      "AGU": "S",      "GGU": "G",
    "UGC": "C",      "CGC": "R",      "AGC": "S",      "GGC": "G",
    "UGA": "Stop",   "CGA": "R",      "AGA": "R",      "GGA": "G",
    "UGG": "W",      "CGG": "R",      "AGG": "R",      "GGG": "G" 
    }

def translate(RNA_str):
    """
Protein Translation
http://rosalind.info/problems/prot/



    The Genetic Codeclick to collapse
    Figure 1. The human hemoglobin molecule consists of 4 polypeptide chains; 伪 subunits are shown in red and 尾 subunits are shown in blue

    Just as nucleic acids are polymers of nucleotides, proteins are chains of smaller molecules called amino acids; 20 amino acids commonly appear in every species. Just as the primary structure of a nucleic acid is given by the order of its nucleotides, the primary structure of a protein is the order of its amino acids. Some proteins are composed of several subchains called polypeptides, while others are formed of a single polypeptide; see Figure 1.

    Proteins power every practical function carried out by the cell, and so presumably, the key to understanding life lies in interpreting the relationship between a chain of amino acids and the function of the protein that this chain of amino acids eventually constructs. Proteomics is the field devoted to the study of proteins.

    How are proteins created? The genetic code, discovered throughout the course of a number of ingenious experiments in the late 1950s, details the translation of an RNA molecule called messenger RNA (mRNA) into amino acids for protein creation. The apparent difficulty in translation is that somehow 4 RNA bases must be translated into a language of 20 amino acids; in order for every possible amino acid to be created, we must translate 3-nucleobase strings (called codons) into amino acids. Note that there are 43=64 possible codons, so that multiple codons may encode the same amino acid. Two special types of codons are the start codon (AUG), which codes for the amino acid methionine always indicates the start of translation, and the three stop codons (UAA, UAG, UGA), which do not code for an amino acid and cause translation to end.

    The notion that protein is always created from RNA, which in turn is always created from DNA, forms the central dogma of molecular biology. Like all dogmas, it does not always hold; however, it offers an excellent approximation of the truth.

    A eukaryotic organelle called a ribosome creates peptides by using a helper molecule called transfer RNA (tRNA). A single tRNA molecule possesses a string of three RNA nucleotides on one end (called an anticodon) and an amino acid at the other end. The ribosome takes an RNA molecule transcribed from DNA (see 鈥淭ranscribing DNA into RNA鈥??), called messenger RNA (mRNA), and examines it one codon at a time. At each step, the tRNA possessing the complementary anticodon bonds to the mRNA at this location, and the amino acid found on the opposite end of the tRNA is added to the growing peptide chain before the remaining part of the tRNA is ejected into the cell, and the ribosome looks for the next tRNA molecule.

    Not every RNA base eventually becomes translated into a protein, and so an interval of RNA (or an interval of DNA translated into RNA) that does code for a protein is of great biological interest; such an interval of DNA or RNA is called a gene. Because protein creation drives cellular processes, genes differentiate organisms and serve as a basis for heredity, or the process by which traits are inherited.

Problem

The 20 commonly occurring amino acids are abbreviated by using 20 letters from the English alphabet (all letters except for B, J, O, U, X, and Z). Protein strings are constructed from these 20 symbols. Henceforth, the term genetic string will incorporate protein strings along with DNA strings and RNA strings.

The RNA codon table dictates the details regarding the encoding of specific codons into the amino acid alphabet.

Given: An RNA string s corresponding to a strand of mRNA (of length at most 10 kbp).

Return: The protein string encoded by s.
Sample Dataset

AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA

Sample Output

MAMAPRTEINSTRING

    >>> translate("AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA")
    'MAMAPRTEINSTRING'
"""
    
    RNA_str = RNA_str.rstrip()
    if len(RNA_str)%3:   # len(RNA_str)%3 != 0
        raise Exception("RNA sequence should have a long of 3*N: %s" % RNA_str)

    result = []
    for i in range(len(RNA_str)//3):
        codon = RNA_str[i*3: i*3+3]
        if codon in ("UAA", "UAG", "UGA"):
            if i != len(RNA_str)//3-1:
                print result
                print codon
                raise Exception("RNA sequence (length: %d) not ended with proper stop codon: %s" % (len(RNA_str), RNA_str))
        else:
            result.append(RNA2protein[codon])

    return "".join(result)

#-------------------------------------------------------------------------------


#===============================================================================
if __name__ == '__main__':
    import doctest
    doctest.testmod()

