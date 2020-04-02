
SEQ = 'ATGCRYKMVBDHWNSatgcrykmvbdhwns'
ASEQ = 'TACGYRMKBVHDWNStacgyrmkbvhdwns'

from math import log
from math import floor
# from string import maketrans
# from string import join
import random
import os
location = {}
import os
import sys
import string
import copy
from importlib import reload


DegenerateMappings = {'R':['A','G'],
                      'Y':['C','T'],
                      'M':['A','C'],
                      'K':['G','T'],
                      'S':['C','G'],
                      'W':['A','T'],
                      'H':['A','C','T'],
                      'B':['C','G','T'],
                      'V':['A','C','G'],
                      'D':['A','G','T'],
                      'N':['A','C','G','T',],
                      'A':['A'],
                      'C':['C'],
                      'G':['G'],
                      'T':['T'],
                      }


    
def degSW(sA,sB, match = 1, MM = -3, gapOpen = -5, gapExtend = -2, returnLocationInfo = None, Verbose = None, DEBUG = None):
    """|degenerateSW|
    seqA,seqB,match = 8, MM = -4,gap =-5,Verbose = None|
    Poor Man's smith waterman||"""
    #Bug Report Fixed 062308
    #for some reason this does not work, it does not find the correct sueqB sequence
    #MB.degSW('AGGTTCGAAATATTCCCCAAAGAAAGCTCA','ATTCGAAATATTACCCCAAAGAAAGCTCAT', match = 1, MM = -3, gapOpen = -5, gapExtend = -2,Verbose = 1)
    #('TTCGAAATATTCCCC-AAAGAAAGCTCA', 'TTCGAAATATTACCCAAAAGAAAGCTCA', 19)
    #   TTCGAAATATTACCCCAAAGAAAGCTCAT 4 C's 3 A's that's OK
    #   TTCGAAATATTACCCAAAAGAAAGCTCA 3C's 4 A's
    #   possibly because at that junction, there are equal paths through that matrix point
    # 
    seqA,seqB = [sA.upper().replace(' ',''),sB.upper().replace(' ','')]

    SeqALen = len(seqA)

    SeqBLen = len(seqB)
    seqA = BaseToNum(seqA)
    seqB = BaseToNum(seqB)

    Rmatrix = [ [[0,-1] for i in range(SeqALen + 1)] for j in range(SeqBLen +1)] ##SeqA is columns SeqB is Rows
    scoreMax = 0
    colMax = 0
    rowMax = 0
    
    for row in range(1,SeqBLen+1,1): ##Build matrix of 1's and 0's
        
        baseB = seqB[row-1]
        for col in range(1,SeqALen+1,1):
            m = MM
            g = 0
            fork = []
            baseA = seqA[col-1]
            if baseB & baseA:
                m = match
            direction = 0
            Inval = Rmatrix[row][col-1][0]
            inDir = Rmatrix[row][col-1][1]
            g = gapOpen
            if inDir == 0 or inDir == 2: g = gapExtend
            val = Inval + m + g
            fork.append([val,direction])
            
            direction = 1
            Inval = Rmatrix[row-1][col-1][0]
            inDir = Rmatrix[row-1][col-1][1]
            g = 0
            val = Inval + m + g +1 #this plus 1 lets diagonal paths be favored
            fork.append([val,direction])
            
            direction = 2
            Inval = Rmatrix[row-1][col][0]
            inDir = Rmatrix[row-1][col][1]
            g = gapOpen
            if inDir == 0 or inDir == 2: g = gapExtend
            val = Inval + m + g
            fork.append([val,direction])
            
            fork.sort(key=lambda x:int(x[0]))
            choice = fork[-1]
            if choice[0] < 0: choice[0] = 0
            Rmatrix[row][col] = choice
            score = choice[0]
            if score > scoreMax:
                scoreMax = score
                colMax = col
                rowMax = row
    if DEBUG: MatrixHumanReadable(seqA, seqB, Rmatrix)
    retSeqA = []
    retSeqB = []

    r = rowMax
    c = colMax
    t = 1
    direction = Rmatrix[r][c][1]
    score = Rmatrix[r][c][0]
    while direction != -1 and score != 0:
        if direction == 0:
            retSeqA = [seqA[c-1]] + retSeqA
            retSeqB = ['-'] + retSeqB
            r = r
            c = c -1
        elif direction == 1:
            retSeqA = [seqA[c-1]] + retSeqA
            retSeqB = [seqB[r-1]] + retSeqB
            r = r -1
            c = c -1
        elif direction == 2:
            retSeqA = ['-'] + retSeqA
            retSeqB = [seqB[r-1]] + retSeqB  #BugFixed 6-23-08: retSeqB = [seqB[r]] + retSeqB
            r = r -1
            c = c
        score = Rmatrix[r][c][0]
        direction = Rmatrix[r][c][1]
    retSeqA = NumToBase(retSeqA)
    retSeqB = NumToBase(retSeqB)
    if Verbose and scoreMax > 6:
        print (retSeqA)
        print (retSeqB)
        print ('MaxScore: %i' %(scoreMax))
        print ('Position is from %i to %i in seqB' %(rowMax, r))
        print ('Position is from %i to %i in seqA' %(colMax, c))
    if returnLocationInfo:
        return colMax, c, rowMax, r, retSeqA, retSeqB, scoreMax
    return retSeqA,retSeqB,scoreMax

def BaseToNum(seq):
    seq = seq.upper()
    BasetoNumDict = {'C': 1,
                'T': 2,
                'Y': 3,
                'A': 4,
                'M': 5,
                'W': 6,
                'H': 7,
                'G': 8,
                'S': 9,
                'K': 10,
                'B': 11,
                'R': 12,
                'V': 13,
                'D': 14,
                'N': 15,
                'I': 15,
                     'X':15}
    if len(seq) ==1:
        return BasetoNumDict[seq]
    else:
        baseLst = []
        for base in seq:
            baseLst.append(BasetoNumDict[base])
        return baseLst

def NumToBase(inLst):
    NumToBaseDict = {1: 'C',
                     2: 'T',
                     3: 'Y',
                     4: 'A',
                     5: 'M',
                     6: 'W',
                     7: 'H',
                     8: 'G',
                     9: 'S',
                     10: 'K',
                     11: 'B',
                     12: 'R',
                     13: 'V',
                     14: 'D',
                     15: 'N',
                     '-':'-'}
    if not('list' in str(type(inLst))):
        return NumToBaseDict[inLst]
    outseq = ''
    for base in inLst:
        outseq = outseq + NumToBaseDict[base]
    return outseq


def MatrixHumanReadable(seqA, seqB, M):
    dirdct = {-1:'|',0:'<',1:'\\',2:'^'}
    print ('_\t\t',)
    for b in seqA:##Columns
        print (b + '\t',)
    print ('')
    r = -1
    for row in M:
        if r ==-1:
            print ('i\t',)
        else:
            print (seqB[r] + '\t',)##Rows
        r +=1
        for point in row:
            score, direction = point
            print (dirdct[direction]+ str(score) + '\t',)
        print ('')
    print ('end of Matrix')

    
def MakeMatchLst(hit, labelA = 'seqA', labelB = 'seqB'):
    matchLst = []
    mm = 0
    align_len=len(hit[0])
    for i in range(align_len):
        if hit[0][i] == hit[1][i]:
            matchLst.append('|')
        else:
            mm+=1
            matchLst.append('x')
    matchLine = ''.join(matchLst)
    hitStr = "%s - %s \n%s - matchline \n%s - %s" %(hit[0], labelA, matchLine, hit[1], labelB)
    three_prime_mm_idxLst=[]
    for i in range(len(matchLst)):
        if matchLst[i]=='x':
            three_prime_mm_idxLst.append(i-align_len)
    three_prime_mm_idxLst
    return mm, hitStr, three_prime_mm_idxLst

def SW(sA,sB, match = 1, MM = -3, gapOpen = -5, gapExtend = -2, labelA = 'A', labelB = 'B', Verbose = None, DEBUG = None, ReturnHit=False, ReturnPercHit=False):
    """|SW|
    seqA,seqB,match = 8, MM = -4,gap =-5,Verbose = None|
    Poor Man's smith waterman||"""
    seqA,seqB = [sA.upper(),sB.upper()]

    SeqALen = len(seqA)

    SeqBLen = len(seqB)


    Rmatrix = [ [[0,-1] for i in range(SeqALen + 1)] for j in range(SeqBLen +1)]
    scoreMax = 0
    colMax = 0
    rowMax = 0
    for row in range(1,SeqBLen+1,1): ##Build matrix of 1's and 0's
        
        baseB = seqB[row-1]
        for col in range(1,SeqALen+1,1):
            m = MM
            g = 0
            fork = []
            baseA = seqA[col-1]
            if baseB == baseA:
                m = match
            direction = 0#Left
            Inval = Rmatrix[row][col-1][0]
            inDir = Rmatrix[row][col-1][1]
            g = gapOpen
            if inDir == 0 or inDir == 2: g = gapExtend
            val = Inval + m + g
            fork.append([val,direction])
            
            direction = 1#diagonal
            Inval = Rmatrix[row-1][col-1][0]
            inDir = Rmatrix[row-1][col-1][1]
            g = 0
            val = Inval + m + g
            fork.append([val,direction])
            
            direction = 2#Up
            Inval = Rmatrix[row-1][col][0]
            inDir = Rmatrix[row-1][col][1]
            g = gapOpen
            if inDir == 0 or inDir == 2: g = gapExtend
            val = Inval + m + g
            fork.append([val,direction])
            
            fork.sort(key=lambda x:int(x[0]))
            choice = fork[-1]
            if choice[0] < 0: choice[0] = 0
            Rmatrix[row][col] = choice
            score = choice[0]
            if score > scoreMax:
                scoreMax = score
                colMax = col
                rowMax = row

    retSeqA = ''
    retSeqB = ''

    r = rowMax
    c = colMax
    t = 1
    direction = Rmatrix[r][c][1]
    score = Rmatrix[r][c][0]
    if DEBUG: MatrixHumanReadable(seqA, seqB, Rmatrix)
    while direction != -1 and score != 0:
        if direction == 0:
            retSeqA = seqA[c-1] + retSeqA
            retSeqB = '-' + retSeqB
            r = r
            c = c -1
        elif direction == 1:
            retSeqA = seqA[c-1] + retSeqA
            retSeqB = seqB[r-1] + retSeqB
            r = r -1
            c = c -1
        elif direction == 2:
            retSeqA = '-' + retSeqA
            retSeqB = seqB[r-1] + retSeqB#old version: retSeqB = seqB[r] + retSeqB
            r = r -1
            c = c
        score = Rmatrix[r][c][0]
        direction = Rmatrix[r][c][1]
    if ReturnHit:
        align_len=len(retSeqA)
        mm, match_str, three_prime_mm_idxLst=MakeMatchLst((retSeqA, retSeqB), labelA = labelA, labelB = labelB)
        three_prime_mm=''
        if three_prime_mm_idxLst:
            three_prime_mm=max(three_prime_mm_idxLst)
        three_prime_overhang=len(sA)-c-align_len
        align_percentHit = 1.0- (mm/align_len if align_len>0 else 0)
        hit={'a_end':colMax,
           'b_end':rowMax,
           'a_start':c,
           'b_start':r,
           'retSeqA':retSeqA,
           'retSeqB':retSeqB,
           'match_str':match_str,
           'num_mismatches':mm,
           'align_len':align_len,
           'align_percentHit':align_percentHit,
           'a_percHit':(align_len-mm)/len(sA)  if len(sA)>0 else 0,
           'b_percHit':(align_len-mm)/len(sB)  if len(sB)>0 else 0,
           'three_prime_mm':three_prime_mm,
           'three_prime_overhang':three_prime_overhang,
           'test':(sA[c:colMax], sB[r:rowMax])
          }
        return hit
    if ReturnPercHit:
        mm, match_str, three_prime_mm_idxLst=MakeMatchLst((retSeqA, retSeqB), labelA = 'A', labelB = 'B')
        align_len=len(retSeqA)
        return [(align_len-mm)/len(sA), (align_len-mm)/len(sB)]
    if Verbose:
        print (retSeqA)
        print (retSeqB)
        print ('MaxScore: %i' %(scoreMax))
    return [retSeqA,retSeqB]
        

def computeFreeEnergy(p_line, m_line, subRoutineDir, Gap = None, gapOpen = -3):
    """computeFreeEnergy
    |Args:
    |Output:
    |Desc:
    |Input:Ex
    |Aux: """
    if not(Gap) and len(p_line) == len(m_line):
        queryHitSeq, subjectHitSeq = [p_line, m_line]
    else:
        queryHitSeq, subjectHitSeq = SW(p_line,m_line, gapOpen = gapOpen, Verbose = None)
    f= os.popen(subRoutineDir + 'energy ' + queryHitSeq + ' ' + subjectHitSeq + ' 1')
    ##In reverting you are correctly antisensing the oppsoite strand. It is not
    ## Real antisensing as the sequence is not flipped it is simply changed in place
    ## GATTC becomes CTAAG
    retData = f.read()
    f.close()
    energy =  float(retData.split('\n') [-1])
    return energy,queryHitSeq, subjectHitSeq
##def computeFreeEnergy(p_line, m_line, subRoutineDir, mm = -3, gapopen = -5, gapextend = -2, Gap = None):
##    #This requires Derisis's energy program in subroutine directory [send this 'c:\pythonCore\']look at Derisi's program oligoarray
##    if not(Gap) and len(p_line) == len(m_line):
##        queryHitSeq, subjectHitSeq = [p_line, m_line]
##    else:
##        queryHitSeq, subjectHitSeq = SW(p_line,m_line, MM = mm, gapOpen = gapopen, gapExtend = gapextend, Verbose = None)
##    f= os.popen(subRoutineDir + 'energy ' + queryHitSeq + ' ' + subjectHitSeq + ' 1')
##    ##In reverting you are correctly antisensing the oppsoite strand. It is not
##    ## Real antisensing as the sequence is not flipped it is simply changed in place
##    ## GATTC becomes CTAAG
##    retData = f.read()
##    f.close()
##    energy =  float(retData.split('\n') [-1])
##    return energy,queryHitSeq, subjectHitSeq

def computeFreeEnergyHack(queryHitSeq, subjectHitSeq, location=None):
##    subRoutineDir = location['subRoutineDir']
##    UD = location['WorkDir']
##    f= os.popen(subRoutineDir + 'energy ' + queryHitSeq + ' ' + subjectHitSeq + ' 1')
##    ##In reverting you are correctly antisensing the oppsoite strand. It is not
##    ## Real antisensing as the sequence is not flipped it is simply changed in place
##    ## GATTC becomes CTAAG
##    retData = f.read()
##    f.close()
##    if 'WARNING:' in retData:
##        energy = 0.0
##    else:
##        energy =  float(retData.split('\n') [-1])
    qLst = queryHitSeq.split('-')
    qSum = 0.0
    for s in qLst:
        qSum += dG(s)

    sLst = subjectHitSeq.split('-')
    sSum = 0.0
    for s in sLst:
        sSum += dG(s)
    sLst = queryHitSeq.split('-')
    energy = (qSum + sSum)/2
    return energy

def CodonForm(seq):
    i = 0
    outseq = ''
    for s in seq:
        outseq = outseq+s
        i+=1
        if i ==3:
            i = 0
            outseq += ' '
    return outseq

def AAtoNA(base):
    aana = {'A':['GCC','GCA','GCG','GCT'],
            'C':['TGT','TGC'],
            'D':['GAT','GAC'],
            'E':['GAA','GAG'],
            'F':['TTT','TTC'],
            'G':['GGC','GGA','GGG','GGT'],
            'H':['CAT','CAC'],
            'I':['ATC','ATA','ATT'],
            'K':['AAA','AAG'],
            'L':['TTA','CTA','TTG','CTG','CTT','CTG'],
            'M':['ATG'],
            'N':['AAT','AAC'],
            'P':['CCC','CCA','CCG','CCT'],
            'Q':['CAA','CAG'],
            'R':['CGT','CGC','AGA','CGA','AGG','CGA'],
            'S':['TCT','TCC','AGT','TCA','AGC','TCA'],
            'T':['ACT','ACC','ACA','ACG'],
            'V':['GTG','GTT','GTC','GTA'],
            'W':['TGG'],
            'Y':['TAT','TAC']}
    if base.upper() in aana:
        return aana[base]
    else:
        return []
    
def NAtoAA(codon):
    DNAtoAA = { 'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L',
             'CTA': 'L', 'CTG': 'L', 'CTN': 'L', 'TGG': 'W',
             'TAA': '*', 'TAG': '*', 'TGA': '*', 'ATG': 'M',
             'TTT': 'F', 'TTC': 'F', 'TAT': 'Y', 'TAC': 'Y',
             'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 
             'TCN': 'S', 'AGT': 'S', 'AGC': 'S', 'CCT': 'P', 
             'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CCN': 'P',
             'TGT': 'C', 'TGC': 'C', 'CAT': 'H', 'CAC': 'H',
             'CAA': 'Q', 'CAG': 'Q', 'AAT': 'N', 'AAC': 'N',
             'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 
             'CGN': 'R', 'AGA': 'R', 'AGG': 'R', 'ATT': 'I', 
             'ATC': 'I', 'ATA': 'I', 'AAA': 'K', 'AAG': 'K',
             'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
             'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
             'ACN': 'T', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V',
             'GTG': 'V', 'GTN': 'V', 'GCT': 'A', 'GCC': 'A',
             'GCA': 'A', 'GCG': 'A', 'GCN': 'A', 'GGT': 'G', 
             'GGC': 'G', 'GGA': 'G', 'GGG': 'G', 'GGN': 'G',
             'TAN': 'X', 'TTN': 'X', 'TGN': 'X', 'CAN': 'X', 
             'ATN': 'X', 'AAN': 'X', 'GAN': 'X', 'AGN': 'X',
             'ANA': 'X', 'ANT': 'X', 'ANG': 'X', 'ANC': 'X', 
             'TNA': 'X', 'TNT': 'X', 'TNG': 'X', 'TNC': 'X', 
             'GNA': 'X', 'GNT': 'X', 'GNG': 'X', 'GNC': 'X', 
             'CNA': 'X', 'CNT': 'X', 'CNG': 'X', 'CNC': 'X', 
             'NAA': 'X', 'NAT': 'X', 'NAG': 'X', 'NAC': 'X', 
             'NTA': 'X', 'NTT': 'X', 'NTG': 'X', 'NTC': 'X', 
             'NGA': 'X', 'NGT': 'X', 'NGG': 'X', 'NGC': 'X', 
             'NCA': 'X', 'NCT': 'X', 'NCG': 'X', 'NCC': 'X', 
             'NNA': 'X', 'NNT': 'X', 'NNG': 'X', 'NNC': 'X', 
             'ANN': 'X', 'TNN': 'X', 'GNN': 'X', 'NNC': 'X', 
             'NAN': 'X', 'NTN': 'X', 'NGN': 'X', 'NCN': 'X', 
             'NNN': 'X','':''
                }
    return DNAtoAA[codon.upper()]

def NAtoProtein(inseq):
    """Translate|
    (inseq)
    |returns protein sequence a s tring
    |uses a dictionary
    |"""
    ##            """ Ter = ["TAA", "TAG", "TGA"]
    ##            Phe = ["TTT", "TTC"]
    ##            Leu = ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG", "CTN"]
    ##            Ser = ["TCT", "TCC", "TCA", "TCG", "TCN", "AGT", "AGC"]
    ##            Tyr = ["TAT", "TAC"]
    ##            Cys = ["TGT", "TGC"]
    ##            Trp = ["TGG"]
    ##            Pro = ["CCT", "CCC", "CCA", "CCG", "CCN"]
    ##            His = ["CAT", "CAC"]
    ##            Gln = ["CAA", "CAG"]
    ##            Arg = ["CGT", "CGC", "CGA", "CGG", "CGN", "AGA", "AGG"]
    ##            Ile = ["ATT", "ATC", "ATA"]
    ##            Met = ["ATG"]
    ##            Thr = ["ACT", "ACC", "ACA", "ACG", "ACN"]
    ##            Asn = ["AAT", "AAC"]
    ##            Lys = ["AAA", "AAG"]
    ##            Val = ["GTT", "GTC", "GTA", "GTG", "GTN"]
    ##            Ala = ["GCT", "GCC", "GCA", "GCG", "GCN"]
    ##            Asp = ["GAT", "GAC"]
    ##            Glu = ["GAA", "GAG"]
    ##            Gly = ["GGT", "GGC", "GGA", "GGG", "GGN"]
    ##            Any = ["TAN", "TTN", "TGN", "CAN", "ATN", "AAN", "GAN", "AGN", \
    ##                           "ANA", "ANT", "ANG", "ANC", \
    ##                           "TNA", "TNT", "TNG", "TNC", \
    ##                           "GNA", "GNT", "GNG", "GNC", \
    ##                           "CNA", "CNT", "CNG", "CNC", \
    ##                           "NAA", "NAT", "NAG", "NAC", \
    ##                           "NTA", "NTT", "NTG", "NTC", \
    ##                           "NGA", "NGT", "NGG", "NGC", \
    ##                           "NCA", "NCT", "NCG", "NCC", \
    ##                           "NNA", "NNT", "NNG", "NNC", \
    ##                           "ANN", "TNN", "GNN", "NNC", \
    ##                           "NAN", "NTN", "NGN", "NCN", \
    ##                           "NNN"]
    ##            Any : "X", Ter: "*", Phe: "F",
    ##            Leu: "L", Ser: "S", Tyr: "Y",
    ##            Cys: "C", Trp: "W", Pro: "P",
    ##            His: "H", Gln: "Q", Arg: "R",
    ##            Ile: "I", Met: "M", Thr: "T",	   
    ##            Asn: "N", Lys: "K", Val: "V",
    ##            Ala: "A", Asp: "D", Glu: "E",
    ##            Gly: "G" """

    DNAtoAA = { 'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L',
             'CTA': 'L', 'CTG': 'L', 'CTN': 'L', 'TGG': 'W',
             'TAA': '*', 'TAG': '*', 'TGA': '*', 'ATG': 'M',
             'TTT': 'F', 'TTC': 'F', 'TAT': 'Y', 'TAC': 'Y',
             'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 
             'TCN': 'S', 'AGT': 'S', 'AGC': 'S', 'CCT': 'P', 
             'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CCN': 'P',
             'TGT': 'C', 'TGC': 'C', 'CAT': 'H', 'CAC': 'H',
             'CAA': 'Q', 'CAG': 'Q', 'AAT': 'N', 'AAC': 'N',
             'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 
             'CGN': 'R', 'AGA': 'R', 'AGG': 'R', 'ATT': 'I', 
             'ATC': 'I', 'ATA': 'I', 'AAA': 'K', 'AAG': 'K',
             'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
             'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
             'ACN': 'T', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V',
             'GTG': 'V', 'GTN': 'V', 'GCT': 'A', 'GCC': 'A',
             'GCA': 'A', 'GCG': 'A', 'GCN': 'A', 'GGT': 'G', 
             'GGC': 'G', 'GGA': 'G', 'GGG': 'G', 'GGN': 'G',
             'TAN': 'X', 'TTN': 'X', 'TGN': 'X', 'CAN': 'X', 
             'ATN': 'X', 'AAN': 'X', 'GAN': 'X', 'AGN': 'X',
             'ANA': 'X', 'ANT': 'X', 'ANG': 'X', 'ANC': 'X', 
             'TNA': 'X', 'TNT': 'X', 'TNG': 'X', 'TNC': 'X', 
             'GNA': 'X', 'GNT': 'X', 'GNG': 'X', 'GNC': 'X', 
             'CNA': 'X', 'CNT': 'X', 'CNG': 'X', 'CNC': 'X', 
             'NAA': 'X', 'NAT': 'X', 'NAG': 'X', 'NAC': 'X', 
             'NTA': 'X', 'NTT': 'X', 'NTG': 'X', 'NTC': 'X', 
             'NGA': 'X', 'NGT': 'X', 'NGG': 'X', 'NGC': 'X', 
             'NCA': 'X', 'NCT': 'X', 'NCG': 'X', 'NCC': 'X', 
             'NNA': 'X', 'NNT': 'X', 'NNG': 'X', 'NNC': 'X', 
             'ANN': 'X', 'TNN': 'X', 'GNN': 'X', 'NNC': 'X', 
             'NAN': 'X', 'NTN': 'X', 'NGN': 'X', 'NCN': 'X', 
             'NNN': 'X'
                }
    proteinSeq = ''
    position = 0
    errors = []
    for c in range(0,len(inseq),3):
        position +=1
        codon = inseq[c:c+3].upper()
        if codon in DNAtoAA:
            AA = DNAtoAA[codon]
            proteinSeq = proteinSeq+AA
        else:
            errors.append([codon,position])
    return proteinSeq, errors
      
def CountN(inseq, N = 'N'):
    """countN(inseq)|||There must be a better way to do this|"""
    count = 0
    for i in inseq.upper():
        if i == N:
            count +=1
    return count

def CountRE(inseq, RE = 'TCTAGA'):
    """countN(inseq)|||There must be a better way to do this|"""
    count = 0
    inseq = inseq.upper()
    RE = RE.upper()
    pos = inseq.find(RE)
    if pos != -1:
        count +=1
    while pos != -1:
        pos = inseq.find(RE, pos + 1)
        if pos != -1:
            count +=1
    return count


def onebaserepeats(seq):
    repeatMax = 0
    for i in range(len(seq)):
        numrepeats = 1
        for j in range(i,len(seq)-1,1):
            if seq[j] == seq[j+1]:
                
                numrepeats +=1
                if numrepeats > repeatMax:
                    repeatMax = numrepeats
                    REPEATPOS = i
                #print seq[j] + '  ' + seq[j+1] + '  ' + str(numrepeats) + '  '+str(repeatMax)
            else:
                #print 'Break'
                numrepeats = 0
                REPEATPOS = i
                break
    return repeatMax

def repeats(seq, singleOut = None):
    """repeats
    |seq
    |integer
    |determines how many 1*(multiply by 2) or two base repeats in seq
    |you should dump probes with return >3
    """
    seq = seq.upper()
    pieceALst = []
    pieceBLst = []
    maxA = 0
    maxB = 0
    previousA = '  '
    previousB = '  '
    maxApiece = '  '
    maxBpiece = '  '
    incr = 1
    for i in range(0,len(seq),2):
        pieceALst.append(seq[i:i+2])
        pieceBLst.append(seq[i+1:i+3])
    #print str(pieceALst)
    #print str(pieceBLst)
    for obj in pieceALst:
        #print previousA + '  ' + obj
        if obj == previousA:
            incr +=1
            #print 'Incr' + previousA + '  ' + obj
            #print str(incr)
            if incr > maxA:
                #print '%i incr  %i maxA %s' %(incr,maxA, previousA)
                maxA = incr
                maxApiece = previousA
        else:
            incr = 1
        previousA = obj
    incr = 1
    #print 'maxA' + str(maxA)
    for obj in pieceBLst:
        if previousB ==obj:
            incr +=1
            if incr > maxB:
                #print '%i incr  %i maxB %s' %(incr,maxB, previousB)
                maxB = incr
                maxBpiece = previousB
        else:
            incr = 1
        previousB = obj
    #print 'maxB' + str(maxB)
    if maxA > maxB:
        maxLen =  maxA
        maxPiece = maxApiece
    else:
        maxLen =  maxB
        maxPiece = maxBpiece

    oneBaseRepeats = float(onebaserepeats(seq))
    if oneBaseRepeats/2 >= maxLen:
        maxLen = oneBaseRepeats + 0.1
    else:
        maxLen = float(maxLen) + 0.2

    return maxLen
        
    

def dH(sSeq, init = 1):
    """dH
    |(seq, init = 1)
    |Output: -154.500
    |Computes Enthalpy of DNA sequence.
    |Input: 'AAAAAACGTAAACGTACCCC', initiation? (1,0)|"""
    sSeq = sSeq.upper()
    for base in sSeq:
        if not(base in 'acgtACGT'):
    ##if ('N' in sSeq) or sSeq == '':
            return 0
    m = len(sSeq)
    sSeq = sSeq.upper()
    dH = 0
    #AMA Initiation from Allawi-SantaLucia97
    if init:
        #AMA initiation from first base
        initH = {'A':2.3,'T':2.3,'C':0.1,'G':0.1}
        FIRST = sSeq[0:1]
        dH = initH[FIRST]
        #AMA initiation from last base
        lastH = {'A':2.3,'T':2.3,'C':0.1,'G':0.1}
        LAST = sSeq[-1:]
        dH = dH + lastH[LAST]
    enthalpy = {'A':{'A':-7.9,'T':-7.2,'C':-8.4,'G':-7.8},#
                'T':{'A':-7.2,'T':-7.9,'C':-8.2,'G':-8.5},#
                'C':{'A':-8.5,'T':-7.8,'C':-8,'G':-10.6}, #
                'G':{'A':-8.2,'T':-8.4,'C':-9.8,'G':-8}}
    for x in range(len(sSeq)-1):
        Na = sSeq[x:x+1]
        Nb = sSeq[x+1:x+2]
        dH = dH + enthalpy[Na][Nb]
    return dH

def dS(sSeq, init = 1):
    """dS
    |(sSeq, init = 1)
    |Output: -422.7999
    |Computes Entropy of DNA sequence.
    |Input: 'AAAAAACGTAAACGTACCCC', initiation? (1,0)|"""
    sSeq = sSeq.upper()
    for base in sSeq:
        if not(base in 'acgtACGT'):
    ##if ('N' in sSeq) or sSeq == '':
            return 0
    m = len(sSeq)
    sSeq = sSeq.upper()
    dS = 0
     #AMA Initiation from Allawi-SantaLucia97
    if init:
        #AMA initiation from first base
        initS = {'A':4.1,'T':4.1,'C':-2.8,'G':-2.8}
        FIRST = sSeq[0:1]
        dS = initS[FIRST]
        #AMA initiation from last base
        lastS = {'A':4.1,'T':4.1,'C':-2.8,'G':-2.8}
        LAST = sSeq[-1:]
        dS = dS + lastS[LAST]
    entropy = {'A':{'A':-22.2,'T':-20.4,'C':-22.4,'G':-21},#
               'T':{'A':-21.3, 'T':-22.2,'C':-22.2,'G':-22.7},#
               'C':{'A':-22.7,'T':-21,'C':-19.9, 'G':-27.2},#
               'G':{'A':-22.2,'T':-22.4,'C':-24.4,'G':-19.9}}
    for x in range(len(sSeq)-1):
        Na = sSeq[x:x+1]
        Nb = sSeq[x+1:x+2]
        dS = dS + entropy[Na][Nb]
    return dS

def dG(sSeq, init = 1, T = 37.0):
    """dG
    |(sSeq, init = 1, T = 37.0)
    |Output: -23.36858
    |Computes Energy of DNA sequence.
    |Input: 'AAAAAACGTAAACGTACCCC', 1/0, T = 37.0C|"""
    sSeq = sSeq.upper()
    if not('N' in sSeq) and sSeq != '':
        Tkelvin = T + 273.15
        #print str(dH(sSeq, init)) + '  S' + str( dS(sSeq,init))
        deltaG = dH(sSeq, init) - dS(sSeq,init)/1000*Tkelvin
    else:
        deltaG = 0
    return deltaG

def HybTm(seq):
    return Tm(seq, SaltNa = 0.333, SaltMag = 0.0, Conc = 1.0e-9)

def PCRTm(seq):
    return Tm(seq, SaltNa = 0.05, SaltMag = 1.0, Conc = 200.0e-9)

def Tm(sSeq, SaltNa = 0.333, SaltMag = 0.0, Conc = 1e-9):
    """Tm
    |(sSeq, SaltNa = 0.333, SaltMag = 0.0, Conc = 1e-9)
    |Output:57.8753
    |Computes Tm of DNA sequence at microarray conditions.
    |Input: 'AAAAAACGTAAACGTACCCC', (Na), (Mg), (Probe)|"""
    T = 310.15
    R = 1.987
    SaltMag = float(SaltMag)
    sSeq = sSeq.upper()
    for base in sSeq:
        if not(base in 'GATCgatc'):
            return 0
    if sSeq != '':
        dHh = dH(sSeq)
        dSs = dS(sSeq)
        salt = SaltNa + 4 * pow((SaltMag/ 1000), 0.5)
        saltCorrection = 16.6 * (log(1.7) + log(salt / (1.0 + 0.7 * salt))) / log(10)
        BetterTm = dHh * 1000.0 / (dSs + R * log(Conc / 4.0)) - 273.15 + saltCorrection
        return BetterTm
    
def Antisense(sSeq):
    """Antisense
    |(sSeq)
    |Output:'GGGGTACGTTTACGTTTTTT'
    |Returns reverse complement of any give sequence.
    |Input: 'AAAAAACGTAAACGTACCCC'|"""
    sComp = sSeq.translate(str.maketrans(SEQ,ASEQ))
    LComp = list(sComp)
    LComp.reverse()
    return ''.join(LComp)

def Complement(sSeq):
    """Complements does not Antisense
    |(sSeq)
    |Output:'GGGGTACGTTTACGTTTTTT'
    |Returns reverse complement of any give sequence.
    |Input: 'AAAAAACGTAAACGTACCCC'|"""
    sComp = sSeq.translate(str.maketrans(SEQ,ASEQ))
    
    return sComp

def StrReverse(string):
    Lstr = list(string)
    Lstr.reverse()
    return ''.join(Lstr)

def perGC(sSeq):
    """perGC
    |(sSeq)
    |Output: 40.0
    |Computes percent GC of DNA sequence.
    |Input: 'AAAAAACGTAAACGTACCCC'|"""
    if len(sSeq) == 0:
        return 0
    gc = 0
    sSeq = sSeq.upper()
    for s in sSeq:
        if s in "GC": gc = gc + 1.0
    return (gc*100/len(sSeq))

def countGC(sSeq):
    """perGC
    |(sSeq)
    |Output: 40.0
    |Computes percent GC of DNA sequence.
    |Input: 'AAAAAACGTAAACGTACCCC'|"""
    if len(sSeq) == 0:
        return 0
    gc = 0
    sSeq = sSeq.upper()
    for s in sSeq:
        if s in "GC": gc = gc + 1.0
    return gc

def SSseqExtract(obj,seq):
    ##This is useless!!!
    row = obj[1]
    col = obj[2]
    val = obj[0]
    seqStart = row - val +1
    seqEnd = row +1
    seq = seq[seqStart:seqEnd]
    return [val,seqStart,seqEnd,seq]


def SecStruct(Seq, MinAnneal = 4):
    """SecStruct
    |(sSeq, iMinAnneal = 4)
    |Output: [['ACGT', 4, 2]] = [[hairpin seq, hairpinlength, looplength]]
    |This is done by dynamic programming-smith-waterman-type of search.
    |Input:'AAAACGTAAACGTAAA'|"""
    Seq = Seq.upper()
    SeqLen = len(Seq)
    ASeq = Antisense(Seq)
    s_axis = list(Seq)
    a_axis = list(ASeq)
    Msize = len(s_axis)##Change this if you actually want to do seq comparison
    Rmatrix = []
    for col in range(Msize): ##Build matrix of 1's and 0's
        res = []
        for row in range(Msize):
            val = 0
            sbase = s_axis[col]
            abase = a_axis[row]
            if abase == sbase:
                val = 1
            res.append(val)
        Rmatrix.append(res)
    for col in range(1,Msize): ##Score 1,s and 0,1 to scores by summing diagonals
        for row in range(1,Msize):
            val = Rmatrix[col][row]
            if val:
                Rmatrix[col][row] = Rmatrix[col-1][row-1] + val #to allow gaps simply sum max of [-1,-1],[0,-1],[-1,0]
    foundArray = []
    Diagonals = []##This will speed things up maybe but it is very gap unfriendly
    for TopDiagonal in range(Msize):
        miniDiag = []
        row = TopDiagonal
        col = 0
        while row < Msize:
            val = Rmatrix[row][col]
            if val >= MinAnneal:
                miniDiag.append([val,row,col])
            row +=1
            col +=1
        Diagonals.append(miniDiag)
    for BottomDiagonal in range(Msize-2,0,-1): #=-1 is to not fall off the list, the other -1 is to avoid duplication
        miniDiag = []
        row = BottomDiagonal
        col = Msize-1
        while row > 0:
            val = Rmatrix[row][col]
            if val > MinAnneal:
                miniDiag.append([val,row,col])
            row -=1
            col -=1
        Diagonals.append(miniDiag)
    HPins = []
    PDs = []
    for Diag in Diagonals:
        Diag.sort(key=lambda x: x[0])
        if len(Diag) > 1:
            if Diag[-1][0] == Diag[-2][0]:
                if Diag[-1][1] > Diag[-2][1]:
                    HPins.append([Diag[-2],Diag[-1]])
                else:
                    HPins.append([Diag[-1],Diag[-2]])
            else:
                PDs.append(Diag[-1])
        if len(Diag) == 1:
            PDs.append(Diag[-1])
    finalHP = []
    for HP in HPins:
        obj = HP[0]
        nextobj = HP[1]
        #[val,seqStart,seqEnd,seq]
        obj = SSseqExtract(obj,Seq)
        nextobj =SSseqExtract(nextobj,Seq)
        finalHP.append([obj[3],obj[0],nextobj[1]-obj[2]])
    for PD in PDs:
        PD = SSseqExtract(PD,Seq)
        if PD[3] == Antisense(PD[3]):
            finalHP.append([PD[3],PD[0],0])           
    return finalHP



def OldSecStruct(sSeq, iMinAnneal = 4):
    """SecStruct
    |(sSeq, iMinAnneal = 4)
    |Output: [['ACGT', 4, 2]] = [[hairpin seq, hairpinlength, looplength]]
    |Creates a list of all hairpin loops within a sequence making no
    distinction about stability. Excludes loops greater than 9 and palindromes.
    |Input:'AAAACGTAAACGTAAA'|"""
    sSeq = sSeq.upper()
    iSeqLen = len(sSeq)
    sASeq = Antisense(sSeq)
    iQStart = -1
    iSStart = 0
    lHP = []
    #Find All hits If you remove iMinAnneal, it will get more sensitive
    while iQStart < iSeqLen - iMinAnneal:
        iQStart +=1
        iQEnd = iSeqLen
        while  iQEnd >= iQStart + iMinAnneal:
            sQuery = sSeq[iQStart:iQEnd]
            sSearchSpace = sASeq
            if sSearchSpace[:iQEnd].find(sQuery)!=-1:
                lHP.append([sQuery,iQStart,iQEnd])
                iQStart += iMinAnneal
                break
            else:
                iQEnd += -1
    #Locate All hits
    loutHP = []
    for hit in lHP:
         seq = hit[0]
         Aseq = Antisense(seq)
         sFind = hit[1]
         iEnd = hit[2]
         hitLen = len(seq)
         if sSeq[sFind:].find(Aseq)!=-1:
            ssFind = sSeq.find(Aseq,sFind)
            loop = ssFind - iEnd
            if loop <0: loop = 0
         else:
            ssFind = sFind
            loop = 0
         if seq.find(Aseq)!=-1:
             loutHP.append([seq,hitLen,loop])
         elif loop != 0:
             loutHP.append([seq,hitLen,loop])
                 
    return loutHP


def StabSecStruct(lHPs, OutType= None,Pal = None):
    """StabSecStruct
    |(lHPs, OutType = 0,1)
    |Output: ['ACGTG', 5, 2] or 50.5
    |Returns either most stable secondary structure from list or Tm of most
    stable hairpin.
    |Input: [sequence or SecStruct Output, 1/0]|"""
    chosenHP = ''
    if str(type(lHPs)).find('str')!=-1:
        lHPs = SecStruct(lHPs)
    if len(lHPs) == 0:
        if OutType:
            return [0,'']
        else:
            return 0
    fdGMin = 0.0
    fdG = 0.0
    MaxHP = None
    
    for HP in lHPs:
        useHP = 1
        if not(Pal):
            if HP[2] == 0:
                useHP = None
        if useHP:
            fdG = dG(HP[0])
            if fdG < fdGMin:
                fdGMin = fdG 
                MaxHP = Tm(HP[0],Conc =1)
                chosenHP = HP
    if OutType:
        if len(chosenHP) ==3:
            hpstring =  str(chosenHP[0]) + '_' + str(chosenHP[1]) + '_' + str(chosenHP[2])
        else:
            hpstring = ''
        return [MaxHP , hpstring]
    else: 
        return MaxHP

def Hairpin(seq):
    """Hairpin
    |(seq)
    |Output: [MaxHP , hpstring]
    |Returns output of SecStruct and of StabSecStruct
    ||"""
    probeSeqSSLst = SecStruct(seq, MinAnneal = 4)
    HPLst = []
    for obj in probeSeqSSLst:
        if obj[2]>0:
            HPLst.append(obj)
    HPobj = StabSecStruct(HPLst, OutType = 1)
    #print str(HPobj)
    return HPobj

def HairpinTm(seq, MinAnneal = 4):
    """HairpinTm
    |(seq)
    |Output: HPTm
    |Returns hairpin Tm
    ||"""
    probeSeqSSLst = SecStruct(seq, MinAnneal = MinAnneal)
    HPLst = []
    for obj in probeSeqSSLst:
        if obj[2]>0:
            HPLst.append(obj)
    HPobj = StabSecStruct(HPLst, OutType = 1)
    if HPobj == None:
        return 0
    else:
        if HPobj[0] ==None:
            return 0
        else:
            return HPobj[0]


def TMPCRHYB(seq, PCR):
    if PCR:
        return PCRTm(seq)
    else:
        return HybTm(seq)

def AdjustProbeTmTo(probeSeq, dTm = 70.0, addToFivePrime = 1, RemoveFromFivePrime = 1, Verbose = 1, PCR = 1):
    """AdjustTmTo
    |probeSeq, dTm = 70.0, addToFivePrime = 1, Verbose = 1
    |Output: newprobeSeq
    |Returns longer or shorter Hairpin
    ||"""
    currentTm = TMPCRHYB(probeSeq, PCR)
    if currentTm > dTm:
        #adjust Down To
        count = 0
        while  currentTm > dTm:
            count +=1
            if RemoveFromFivePrime or (Tm(probeSeq[:-1]) > Tm(probeSeq[1:])):#remove the base that has the greatest disruptive ability
                probeSeq = probeSeq[1:] #take from 5' end
            else:
                probeSeq = probeSeq[:-1] #take from 3' end
            currentTm = TMPCRHYB(probeSeq, PCR)
            
            if Verbose:
                print ('%s removed %i bases Tm is %.1f' %(probeSeq, count, currentTm))
        return probeSeq
    else:
        #adjust UP to
        basesToAdd = ['g','c','t','a']
        hpTm = HairpinTm(probeSeq)
        count = 0
        
        while  currentTm < dTm:
            count +=1
            testLst = []
            for base in basesToAdd:
                if addToFivePrime:
                    newprobe = base+probeSeq
                else:
                    newprobe = probeSeq + base
                deltaHP = abs(hpTm - HairpinTm(newprobe, MinAnneal = 3) )
                newTm = TMPCRHYB(newprobe, PCR)
                deltaTm = newTm - currentTm
                if base in ['g','c']:
                    gcscore = random.random()
                else:
                    gcscore = 0
                score = deltaTm - deltaHP + gcscore
                testLst.append([newprobe, deltaHP, deltaTm, score])
            testLst.sort(key=lambda x:x[-1])
            #print str(testLst)
            probeSeq = testLst[-1][0]
            hpdelta = testLst[-1][1]
            currentTm = TMPCRHYB(newprobe, PCR)
            if Verbose:
                print ('%s added %i bases Tm is %.1f  hpdelta : %.2f, score: %.2f' %(probeSeq, count, currentTm, hpdelta, score))
        return probeSeq

        


def CharacterizeProbe(probeObj, FilterParams = None, returnLabel = None):
    """CharacterizeProbe
    |(probeObj, FilterParams = None, returnLabel = None, Verbose = 1)
    |probe object needs to be of the form [label, seq, len, anything else, , , ]
    |
    |"""
    
    if returnLabel:
        return ['ProbeGood','HP','PD','repeatNum','GC','HPstruct','PDstruct','failreason']
    Limits = {'MaxHPTm':55,
            'MaxPalindrome':60,
            'MaxRepeats':6,
            'MaxGC':65.0,
            'MinGC':35.0}
    if FilterParams:
        Limits.update(FilterParams)
    
    
    MaxHPTm = float(Limits['MaxHPTm'])
    MaxPalindrome = float(Limits['MaxPalindrome'])
    MaxRepeats = int(Limits['MaxRepeats'])
    MaxGC = float(Limits['MaxGC'])
    MinGC = float(Limits['MinGC'])
    outLst = []

    failreason = 'good'
    ProbeGood = 1

    probeSeq = probeObj[1]
    for base in probeSeq:
        if not(base in 'GATCgatc'):
            ProbeGood = 0
            failreason = 'badBase'
            break
    if failreason == 'badBase':
        probeObj.extend([0,0,0,0,0,0,0,failreason])
        return probeObj
    
    probeSeqSSLst = SecStruct(probeSeq, MinAnneal = 4)
    HPLst = []
    PDLst = []
    for obj in probeSeqSSLst:
        if obj[2]>0:
            HPLst.append(obj)
        else:
            PDLst.append(obj)
    PDobj = StabSecStruct(PDLst,Pal = 1,OutType = 1)
    PD = PDobj[0]
    PDstruct = PDobj[1]
    HPobj = StabSecStruct(HPLst, OutType = 1)
    HP = HPobj[0]
    HPstruct = HPobj[1]
    if PD > MaxPalindrome:
        ProbeGood = 0
        failreason = 'PD > ' + str(MaxPalindrome)
    if HP > MaxHPTm:
        ProbeGood = 0
        failreason = 'HP > ' + str(MaxHPTm)
    GC = perGC(probeSeq)
    if GC > MaxGC:
        ProbeGood = 0
        failreason = '%GC > ' + str(MaxGC)
    if GC < MinGC:
        ProbeGood = 0
        failreason = '%GC < ' + str(MinGC)
    repeatNum = repeats(probeSeq)
    
    repeatType = int((repeatNum - floor(repeatNum))*10)
    if  repeatNum > MaxRepeats:
        ProbeGood = 0
        failreason = 'repeats > ' + str(MaxRepeats) + '-' + str(repeatType)
    probeObj = probeObj[:7] ##prevents retarded Concatenations
    probeObj.extend([ProbeGood,HP,PD,repeatNum,GC,HPstruct,PDstruct,failreason])
    
    return probeObj


    

def UnitTest():
    seq = 'GATCGGGGGGGgatgctagctgatcgCCCCCCCCCCC'
#     print 'AntiSense: ' + seq + ' -> ' + Antisense(seq)
#     print '%s(%s): %s' %('Complement',seq, Complement(seq))
#     print '%s(%s): %s' %('StrReverse',seq, StrReverse(seq))
#     print '%s(%s): %s' %('perGC',seq, str(perGC(seq)))
#     print '%s(%s): %s' %('countGC',seq, str(countGC(seq)))
#     print '%s(%s): %s' %('StrReverse',seq, StrReverse(seq))
#     
#     print 'SecStruct: ' + str(SecStruct(seq))
#     print 'StabSecStruct: ' + str(StabSecStruct(seq, OutType = 0))
#     print 'StabSecStruct: ' + str(StabSecStruct(seq, OutType = 1))
#     print 'perGC: ' +str(perGC(seq))
#     print 'SW:Smith Wwaterman'
#     out = SW('GATCGTGTGTGT','GATCGTGTGTGT', match = 1, MM = -3, gapOpen = -5, gapExtend = -2,Verbose = 1)
#     print 'SW_Output is:' + str(out)
#     print 'degSW:Degenerate Smith waterman'
#     out = degSW('GATCYTGTGTGT','GATCGTGTGTGT', match = 1, MM = -3, gapOpen = -5, gapExtend = -2,Verbose = 1)
#     print 'degSW_Output is:' + str(out)
#     
#     print 'BaseToNum(): ' + str(BaseToNum('GATTGTGT'))
#     print 'NumToBase(): ' + str(NumToBase([1,2,3,4]))
#     print 'CodonForm(seq): ' + CodonForm('GATGTGTCACTGTATG')
#     print 'AAtoNA(base): ' + str(AAtoNA('Y'))
#     print 'NAtoAA(codon): ' + str(NAtoAA('gat'))
#     print 'NAtoProtein(inseq): ' + str(NAtoProtein('GATGTGTCACTGTATG'))
#     print 'CountN(gatcgtagctgaNNnnnnnNnnnnnngatrcgtag): ' + str(CountN('gatcgtagctgaNNnnnnnNnnnnnngatrcgtag'))
#     print 'CountRE(GATCGGGGGGGgatgtctagatgatcgtctagaCCCCCCCCCCC, RE = TCTAGA): ' + str(CountRE('GATCGGGGGGGgatgtctagatgatcgtctagaCCCCCCCCCCC', RE = 'TCTAGA'))
#     print 'onebaserepeats(ggggatgctagtcgggggggg): ' + str(onebaserepeats('ggggatgctagtcgggggggg'))
#     print 'repeats(ggggatgctagtctatatatatatatatag): ' + str(repeats('ggggatgctagtctatatatatatatata'))
#     seq = 'gatcgtagctagtcgatcggctagtcgat'
#     print 'dH(%s, init = 1): %s' %(seq, dH(seq, init = 1))
#     print 'dS(%s, init = 1): %s' %(seq, dS(seq, init = 1))
#     print 'dG(%s, init = 1): %s' %(seq, dG(seq, init = 1))
#     print 'Tm(%s, , SaltNa = 0.333, SaltMag = 0.0, Conc = 1e-9): %s' %(seq, Tm(seq,SaltNa = 0.333, SaltMag = 0.0, Conc = 1e-9))
#     ##Probe Filter
    seqA = 'gatgctagtcgatcgtagctgatcgatgcta'
    seqB = 'gatcgctgcgctgctcgctgctcgtgctcgctgctcgccccccccccgggg'
    seqC = 'gatcgatcgtagtcgtgatatatcgtagt'
    probeDB = [['probeA', seqA,len(seqA)],
               ['probeB', seqB,len(seqB)],
               ['probeC', seqC,len(seqC)]]
    
    FilterParams = {'MaxHPTm':55,
            'MaxPalindrome':60,
            'MaxRepeats':6,
            'MaxGC':65.0,
            'MinGC':35.0}
    LabelLine = CharacterizeProbe([], FilterParams = {}, returnLabel = 1)
    data = [CharacterizeProbe(probeObj, FilterParams = FilterParams, returnLabel = None) for probeObj in probeDB]
    data = [LabelLine] + data
    print (str(data))
    return data

class logger:
    """logger is called in order to capture all print statements via the
    sysstdout function and redirect them to a log file, You simply
    Hijack the sys.stdout function and replace its screenwriting object with
    a something that writes to both screen and errLogFile
    This is how it needs to be called:
    screenWriter = sys.stdout
    logObj = logger(errFP,screenWriter,Verbose =1)
    At the end, please return sys its function like so:
    sys.stdout = loggerObj.uncouple()
    """
    def __init__(self,FP,screenWriter,Verbose = 1):
        self.fp = FP
        self.verbose = Verbose
        self.screenWriter = screenWriter
    def write(self,s):
        if self.verbose:
            self.screenWriter.write(str(s))
        self.fp.write(str(s))
    def uncouple(self):
        #Call this at the end like so
        #sys.stdout = logObj.uncouple()
        self.fp.close()
        return self.screenWriter
    
if None:
    UD = os.getcwd()
    errFP = open(UD + 'MolceularBiologyModule_log.txt','w')
    screenWriter = sys.stdout
    logObj = logger(errFP,screenWriter,Verbose =1)
    sys.stdout = logObj #this diverts the print statement's stream
    ##ERROR This writes error messages to a file if execution fails
    fsock = open(UD + 'MolceularBiologyModule_error.log', 'w')
    sys.stderr = fsock

    ##This is the part where the unit test is actually called
    d = UnitTest()
    ret = raw_input('It appears that the MoleculaBiology Module Unit test executed OK (any key):')


    sys.stdout = logObj.uncouple()#Returns the print stream to normal so you don't keep writing to the disk
    
    


#################### DEPRECATION LINE ##########################
#################### DEPRECATION LINE ##########################

def FilterProbes(probeLst, FilterParams, InLabelLine = None, returnData = None, Verbose = 1):
    """FilterProbes
    |probeLst, FilterParams, InLabelLine = None, returnData = None, Verbose = 1)
    |Output:[??,??,??,??,'ProbeGood','HP','PD','repeatNum','GC','HPstruct','PDstruct','failreason']
    |Filters through user-settable paramters, a set of probes
    |Input: list of probes in standard format"""
    Limits = {'MaxHPTm':55,
            'MaxPalindrome':60,
            'MaxRepeats':6,
            'MaxGC':65.0,
            'MinGC':35.0}
    Limits.update(FilterParams)

    LabelList = ['ProbeGood','HP','PD','repeatNum','GC','HPstruct','PDstruct','failreason']
    if InLabelLine:
        LabelList = InLabelLine
    else:
        ##Just in case you would ever want to split this up
        for i in len(probeLst[0]):
            LabelList = ['??'] + LabelList
               
    MaxHPTm = float(Limits['MaxHPTm'])
    MaxPalindrome = float(Limits['MaxPalindrome'])
    MaxRepeats = int(Limits['MaxRepeats'])
    MaxGC = float(Limits['MaxGC'])
    MinGC = float(Limits['MinGC'])
    outLst = []
    numINProbes = len(probeLst)
    numOutProbes = 0
    i = 0
    for j in range(numINProbes):
        failreason = 'good'
        ProbeGood = 1
        probeObj = probeLst[j][:]
        probeSeq = probeObj[1]
        for base in probeSeq:
            if not(base in 'GATCgatc'):
                ProbeGood = 0
                failreason = 'badBase'
                break
        if failreason == 'badBase':
            if returnData:
                probeObj.extend([ProbeGood,0,0,0,0,'','',failreason])
                probeLst[j] = probeObj
            continue
        probeSeqSSLst = SecStruct(probeSeq, MinAnneal = 4)
        HPLst = []
        PDLst = []
        for obj in probeSeqSSLst:
            if obj[2]>0:
                HPLst.append(obj)
            else:
                PDLst.append(obj)
        PDobj = StabSecStruct(PDLst,Pal = 1,OutType = 1)
        PD = PDobj[0]
        PDstruct = PDobj[1]
        HPobj = StabSecStruct(HPLst, OutType = 1)
        HP = HPobj[0]
        HPstruct = HPobj[1]
        if PD > MaxPalindrome:
            ProbeGood = 0
            failreason = 'PD > ' + str(MaxPalindrome)
        if HP > MaxHPTm:
            ProbeGood = 0
            failreason = 'HP > ' + str(MaxHPTm)
        GC = perGC(probeSeq)
        if GC > MaxGC:
            ProbeGood = 0
            failreason = '%GC > ' + str(MaxGC)
        if GC < MinGC:
            ProbeGood = 0
            failreason = '%GC < ' + str(MinGC)
        repeatNum = repeats(probeSeq)
        
        repeatType = int((repeatNum - floor(repeatNum))*10)
        if  repeatNum > MaxRepeats:
            ProbeGood = 0
            failreason = 'repeats > ' + str(MaxRepeats) + '-' + str(repeatType)
        probeObj.extend([ProbeGood,HP,PD,repeatNum,GC,HPstruct,PDstruct,failreason])
        
        if returnData:
            probeLst[j] = probeObj
        else:
            if ProbeGood:
                outLst.append(probeObj)
        if ProbeGood:
            numOutProbes += 1
        
        if Verbose and (i==0 or i==500):
            print ('FilterProbes: Processed %i probes of %i. Chose %i probes so far' %(j,numINProbes,numOutProbes))
            i = 0
        i += 1
    if Verbose:
        print ('FilterProbes: got %i probes and accepted %i probes' %(numINProbes,numOutProbes))
        print ('with the following Parameters MaxHPTm = %i , MaxPalindrome = %i , MaxRepeats = %i , MaxGC = %g , MinGC = %g' %(MaxHPTm, MaxPalindrome, MaxRepeats, MaxGC, MinGC))
    if returnData:
        if InLabelLine:
            probeLst =  [LabelList] + probeLst
        return probeLst
    else:
        return outLst
