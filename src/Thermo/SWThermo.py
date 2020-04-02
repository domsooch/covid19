'''
Created on Nov 9, 2010

@author: dsuciu
'''
import os, sys, math
from importlib import reload
sys.path.append('../')
sys.path.append('./')
#from string import maketrans
from . import TwodArray as TD
from . import MolBio as MB
Verbose = None
reload(TD)
#For antisense
SEQ = 'ATGCRYKMVBDHWNSatgcrykmvbdhwns-'
ASEQ = 'TACGYRMKBVHDWNStacgyrmkbvhdwns-'
makeTran = str.maketrans(SEQ,ASEQ)


# An alternative to all this is here: http://biopython.org/DIST/docs/api/Bio.SeqUtils.MeltingTemp-module.html

dH = [['', 'A/T', 'T/A', 'C/G', 'G/C', 'T/G', 'G/T', 'A/G', 'G/A', 'A/C', 'C/A', 'C/T', 'T/C', 'A/A', 'T/T', 'C/C', 'G/G', ],
        ['A/T', '-7.9', '-7.2', '-8.4', '-7.8', '-2.5', '1', '-0.6', '-0.7', '2.3', '5.3', '0.7', '-1.2', '1.2', '-2.7', '0', '-3.1', ],
        ['T/A', '-7.2', '-7.9', '-8.2', '-8.5', '-1.3', '-0.1', '0.7', '3', '3.4', '7.6', '1.2', '1', '4.7', '0.2', '6.1', '1.6', ],
        ['C/G', '-8.5', '-7.8', '-8', '-10.6', '-2.8', '-4.1', '-0.7', '-4', '1.9', '0.6', '-0.8', '-1.5', '-0.9', '-5', '-1.5', '-4.9', ],
        ['G/C', '-8.2', '-8.4', '-9.8', '-8', '-4.4', '3.3', '-0.6', '0.5', '5.2', '-0.7', '2.3', '5.2', '-2.9', '-2.2', '3.6', '-6', ],
        ['T/G', '-0.1', '1', '3.3', '-4.1', '5.8', '-1.4', '', '', '', '', '', '', '', '', '', '', ],
        ['G/T', '-1.3', '-2.5', '-4.4', '-2.8', '4.1', '5.8', '', '', '', '', '', '', '', '', '', '', ],
        ['A/G', '3', '-0.7', '0.5', '-4', '', '', '', '', '', '', '', '', '', '', '', '', ],
        ['G/A', '0.7', '-0.6', '-0.6', '-0.7', '', '', '', '', '', '', '', '', '', '', '', '', ],
        ['A/C', '7.6', '5.3', '-0.7', '0.6', '', '', '', '', '', '', '', '', '', '', '', '', ],
        ['C/A', '3.4', '2.3', '5.2', '1.9', '', '', '', '', '', '', '', '', '', '', '', '', ],
        ['C/T', '1', '-1.2', '5.2', '-1.5', '', '', '', '', '', '', '', '', '', '', '', '', ],
        ['T/C', '1.2', '0.7', '2.3', '-0.8', '', '', '', '', '', '', '', '', '', '', '', '', ],
        ['A/A', '4.7', '1.2', '-2.9', '-0.9', '', '', '', '', '', '', '', '', '', '', '', '', ],
        ['T/T', '0.2', '-2.7', '-2.2', '-5', '', '', '', '', '', '', '', '', '', '', '', '', ],
        ['C/C', '6.1', '0', '3.6', '-1.5', '', '', '', '', '', '', '', '', '', '', '', '', ],
        ['G/G', '1.6', '-3.1', '-6', '-4.9', '', '', '', '', '', '', '', '', '', '', '', '', ],
        ]
dS = [['', 'A/T', 'T/A', 'C/G', 'G/C', 'T/G', 'G/T', 'A/G', 'G/A', 'A/C', 'C/A', 'C/T', 'T/C', 'A/A', 'T/T', 'C/C', 'G/G', ],
        ['A/T', '-22.2', '-20.4', '-22.4', '-21', '-8.3', '0.9', '-2.3', '-2.3', '4.6', '14.6', '0.2', '-6.2', '1.7', '-10.8', '-4.4', '-9.5', ],
        ['T/A', '-21.3', '-22.2', '-22.2', '-22.7', '-5.3', '-1.7', '0.7', '7.4', '8', '20.2', '0.7', '0.7', '12.9', '-1.5', '16.4', '3.6', ],
        ['C/G', '-22.7', '-21', '-19.9', '-27.2', '-8', '-11.7', '-2.3', '-13.2', '3.7', '-0.6', '-4.5', '-6.1', '-4.2', '-15.8', '-7.2', '-15.3', ],
        ['G/C', '-22.2', '-22.4', '-24.4', '-19.9', '-12.3', '10.4', '-1', '3.2', '14.2', '-3.8', '5.4', '13.5', '-9.8', '-8.4', '8.9', '-15.8', ],
        ['T/G', '-1.7', '0.9', '10.4', '-11.7', '16.3', '-6.2', '', '', '', '', '', '', '', '', '', '', ],
        ['G/T', '-5.3', '-8.3', '-12.3', '-8', '9.5', '16.3', '', '', '', '', '', '', '', '', '', '', ],
        ['A/G', '7.4', '-2.3', '3.2', '-13.2', '', '', '', '', '', '', '', '', '', '', '', '', ],
        ['G/A', '0.7', '-2.3', '-1', '-2.3', '', '', '', '', '', '', '', '', '', '', '', '', ],
        ['A/C', '20.2', '14.6', '-3.8', '-0.6', '', '', '', '', '', '', '', '', '', '', '', '', ],
        ['C/A', '8', '4.6', '14.2', '3.7', '', '', '', '', '', '', '', '', '', '', '', '', ],
        ['C/T', '0.7', '-6.2', '13.5', '-6.1', '', '', '', '', '', '', '', '', '', '', '', '', ],
        ['T/C', '0.7', '0.2', '5.4', '-4.5', '', '', '', '', '', '', '', '', '', '', '', '', ],
        ['A/A', '12.9', '1.7', '-9.8', '-4.2', '', '', '', '', '', '', '', '', '', '', '', '', ],
        ['T/T', '-1.5', '-10.8', '-8.4', '-15.8', '', '', '', '', '', '', '', '', '', '', '', '', ],
        ['C/C', '16.4', '-4.4', '8.9', '-7.2', '', '', '', '', '', '', '', '', '', '', '', '', ],
        ['G/G', '3.6', '-9.5', '-15.8', '-15.3', '', '', '', '', '', '', '', '', '', '', '', '', ],
        ]




def Reverse(seq):
    sLst = list(seq)
    sLst.reverse()
    return ''.join(sLst)

def Antisense(sSeq):
    """Antisense
    |(sSeq)
    |Output:'GGGGTACGTTTACGTTTTTT'
    |Returns reverse complement of any give sequence.
    |Input: 'AAAAAACGTAAACGTACCCC'|"""
    sComp = sSeq.translate(makeTran)
    return Reverse(sComp)

class ThermoVals:
    def __init__(self, dGTm = 37.1, DNAConc = 0.000001, SaltNa = 0.1, SaltMg = 0.0):
        self.dGTm = dGTm
        self.R = 1.987
        
        self.Hyb_saltCorrection = self.setSalt(SaltNa = 0.333, SaltMg = 0.0)
        self.PCR_saltCorrection = self.setSalt(SaltNa = 0.05, SaltMg = 1.0)
        self.DNAConc = DNAConc
        self.dH = TD.TwoD_Array(inData=dH)
        self.dS = TD.TwoD_Array(inData=dS)
        self.dG = TD.TwoD_Array()
        self.keys = self.dH.keys
        #self.CalculatedG(dGTm = self.dGTm)
#    def CalculatedGArray(self, dGTm = 37.1):
#        ##deprecated
#        TmKelvin = 273.15+dGTm
#        dHPairDB = self.dH.PairDB()
#        for pair in dHPairDB:
#            firstKey, secondKey, dHval = pair
#            dSval = self.dS.getVal(firstKey, secondKey)
#            
#            dGval = dHval-dSval/1000.0*TmKelvin
#            self.dG.setVal(firstKey, secondKey, dGval)
    def setSalt(self, SaltNa = 0.1, SaltMg = 0.0):
        self.salt = SaltNa + 4 * math.pow((SaltMg / 1000), 0.5)
        saltCorrection = 16.6 * (math.log(1.7) + math.log(self.salt / (1.0 + 0.7 * self.salt))) / math.log(10)
        return saltCorrection
    def initdH(self, seq):
        initH = {'A':2.3,'T':2.3,'C':0.1,'G':0.1}
        FIRST = seq[0].upper()
        dH = initH[FIRST]
        #AMA initiation from last base
        lastH = {'A':2.3,'T':2.3,'C':0.1,'G':0.1}
        LAST = seq[-1].upper()
        dH = dH + lastH[LAST]
        return dH
    def initdS(self, seq):
        initS = {'A':4.1,'T':4.1,'C':-2.8,'G':-2.8}
        FIRST = seq[0].upper()
        dS = initS[FIRST]
        #AMA initiation from last base
        lastS = {'A':4.1,'T':4.1,'C':-2.8,'G':-2.8}
        LAST = seq[-1].upper()
        dS = dS + lastS[LAST]
        return dS
    def calculatedG(self, dH, dS, dGTm = 37.1):
        TmKelvin = 273.15+dGTm
        dGval = dH-dS/1000.0*TmKelvin
        return dGval
    def Tm(self, seqA, seqB=None, DNAConc = None, saltCorrection = None):
        """Expects seqA is Sense SeqB is Anti-sense"""
        if saltCorrection == None:
            saltCorrection = self.Hyb_saltCorrection
        if DNAConc == None:
            DNAConc = self.DNAConc
        if seqB == None:
            seqB = seqA.translate(makeTran)
        if not(seqB): return 0.0## Fixes issue of sending empty string
        if len(seqA) != len(seqB):
            raise Exception('Two seqs sent to Tm are of different sizes. No Dice!')
        if Verbose:
            print('Tm for:\n%s\n%s' %(seqA, seqB))
        bpLst = []
        for i in range(len(seqA)):
            bp = '%s/%s'%(seqA[i].upper(), seqB[i].upper())
            bpLst.append(bp)
            
        dH = self.initdH(seqA)
        dS = self.initdS(seqA)
        for b in range(len(bpLst)-1):
            firstKey = bpLst[b]
            secondKey = bpLst[b+1]
            dH += float(self.dH.getVal(firstKey, secondKey))
            dS += float(self.dS.getVal(firstKey, secondKey))
            #if Verbose: print '%s - %s : dH= %.2f dS= %.2f' %(firstKey, secondKey, dH, dS)
        dG = self.calculatedG(dH, dS, dGTm = 37.1)
        if dG > 0.0:
            meltTemp = 0.0
        else:
            meltTemp = dH * 1000.0 / (dS + self.R * math.log(DNAConc / 4.0)) - 273.15 + saltCorrection
        if Verbose:
            print('Tm: %.2f\n\n' %meltTemp)
        return meltTemp
    def PureSeqCompare(self, seqA, seqB,labelA = 'A', labelB = 'B'):
        hit = MB.SW(seqA,seqB, match = 1, MM = -1, gapOpen = -2, gapExtend = -1,Verbose = Verbose, DEBUG = 0, ReturnHit=False)
        #You cannot allow gaps into this formula it does not support it yet.
        if Verbose:
            sA = hit[0]
            sB = hit[1]
            matchLst = []
            for i in range(len(sA)):
                if sA[i]==hit[1][i]:
                    matchLst.append('|')
                else:
                    matchLst.append('x')
            print('PureSeqComparison:\t%s\n%s\n%s\n%s\n' %(seqA, sA, ''.join(matchLst), sB))
        return hit
    def HybTm(self, seqA, seqB=None,  DNAConc = 0.000000001):
        if Verbose:
            print('Hyb',)
        return self.Tm(seqA, seqB, DNAConc = DNAConc, saltCorrection =self.Hyb_saltCorrection)
    def PCRTm(self, seqA, seqB=None, DNAConc = 0.00000002):
        if Verbose:
            print('PCR',)
        return self.Tm(seqA, seqB, DNAConc = DNAConc, saltCorrection =self.PCR_saltCorrection)
#    def TmDistance(self, seqA, seqB, DNAConc = None):
#        Tm = self.HybTm(seqA)
#        dTm = self.HybTm(seqA, seqB = seqB)
#        if Verbose:
#            print '%s has Tm %.2f\n when bound to \n%s, it has a Tm of %.2f' %(seqA, Tm, seqB, dTm)
#        return Tm, dTm
    def MakeMatchLst(self, hit, labelA = 'seqA', labelB = 'seqB'):
        #mm, hitStr = self.MakeMatchLst(hit)
        matchLst = []
        mm = 0
        for i in range(len(hit[0])):
            if hit[0][i] == hit[1][i]:
                matchLst.append('|')
            else:
                mm+=1
                matchLst.append('x')
        matchLine = ''.join(matchLst)
        hitStr = "%s - %s \n%s - matchline \n%s - %s" %(hit[0], labelA, matchLine, hit[1], labelB)
        return mm, hitStr
    def SWTm_WillItHyb(self, seqA, seqB, DNAConc = None):
        #Do two sequences in solution interact with each other
        #SeqCompare: SeqA-Sense to SeqB-Reverse.Complement
        #ThermoCompare: SeqA-Hit[0] to SeqB-Hit[1].Complement
        
        hit = self.PureSeqCompare(seqA, Antisense(seqB))
        #ExampleOutput:
        #['GAATTC', 'GAATTC']
        TmA =  self.HybTm(seqA, DNAConc = DNAConc)
        TmB =  self.HybTm(seqB, DNAConc = DNAConc)
        
        if hit[0]:
            sA = hit[0]
            #GAATTC
            #gaattc - ANTISENSE
            #CTTAAG - TRANSLATE
            sB = hit[1].translate(makeTran)
            TmAB = self.HybTm(sA, seqB = sB, DNAConc = DNAConc)
            deltaTm = max(abs(TmA-TmAB), abs(TmB-TmAB))
            return TmA, TmB, TmAB, deltaTm, hit
        else:
            return TmA, 0.0, TmA, 1000, hit
    def SWTm_WillItHyb_json(self, seqA, seqB, DNAConc = None):
        #Do two sequences in solution interact with each other
        #SeqCompare: SeqA-Sense to SeqB-Reverse.Complement
        #ThermoCompare: SeqA-Hit[0] to SeqB-Hit[1].Complement

        hit = self.PureSeqCompare(seqA, Antisense(seqB))
        #ExampleOutput:
        #['GAATTC', 'GAATTC']
        TmA =  self.HybTm(seqA, DNAConc = DNAConc)
        TmB =  self.HybTm(seqB, DNAConc = DNAConc)

        if hit[0]:
            sA = hit[0]
            #GAATTC
            #gaattc - ANTISENSE
            #CTTAAG - TRANSLATE
            sB = hit[1].translate(makeTran)
            TmAB = self.HybTm(sA, seqB = sB, DNAConc = DNAConc)
            deltaTm = max(abs(TmA-TmAB), abs(TmB-TmAB))
            return {
                'TmAB':TmAB,
                'deltaTm':deltaTm,
                'match_str':self.MakeMatchLst(hit, labelA = 'A', labelB = 'B')
            }
            #return TmA, TmB, TmAB, deltaTm, hit
        else:
            return {
                'TmAB':0,
                'deltaTm':1000,
                'match_str':self.MakeMatchLst(hit, labelA = 'A', labelB = 'B')
            }
            #return TmA, 0.0, TmA, 1000, hit
    def SWTm_WillItPCR_primer(self, primer, amplicon, DNAConc=None):
        #Will half the PCR reaction work?
        antisense_amplicon=False
        #hit = self.PureSeqCompare(amplicon, Antisense(primer))
        hit=MB.SW(primer, amplicon, match = 1, MM = -1, gapOpen = -2, gapExtend = -1, labelA = 'Pri', labelB = 'Amp', Verbose = Verbose, DEBUG = 0, ReturnHit=True)
        primerTm =  self.PCRTm(primer, DNAConc = DNAConc)
        # num_mismatches, match_str = self.MakeMatchLst(hit, labelA = 'seq', labelB = 'complement')
        align_score=(float(hit['align_len'])-hit['num_mismatches'])/len(primer)
        
        if align_score>0.6:
            sA = hit['retSeqA']
            sB = hit['retSeqB'].translate(makeTran)
            TmAB = self.PCRTm(sA, seqB=sB, DNAConc=DNAConc)
            deltaTm = primerTm-TmAB
        else:
            deltaTm = primerTm
            TmAB=align_score
        return {
            'primerTm':primerTm,
            'TmAB':TmAB,
            'align_score':align_score,
            'deltaTm':deltaTm,
            'match_str':hit['match_str'],
            'mis_matches':hit['num_mismatches'],
            'hit':hit
        }
    def WillItPCR(self, amplicon, f_primer, r_primer):
        for_wip = self.SWTm_WillItPCR_primer(f_primer, amplicon)
        rev_wip = self.SWTm_WillItPCR_primer(r_primer, MB.Antisense(amplicon))
        ret={}
        if for_wip['TmAB']>0:
            ret['for_Tm']=for_wip['TmAB']
            ret['for_deltaTm']=for_wip['deltaTm']
            ret['for_align_score']=for_wip['align_score']
            ret['for_match']=for_wip['match_str']
            ret['for_mis_matches']=for_wip['mis_matches']
            ret['for_three_prime_mm']=for_wip['hit']['three_prime_mm']
            ret['for_three_prime_overhang']=for_wip['hit']['three_prime_overhang']
        else:
            ret['for_align_score']=0
        if rev_wip['TmAB']>0:
            ret['rev_Tm']=rev_wip['TmAB']
            ret['rev_deltaTm']=rev_wip['deltaTm']
            ret['rev_align_score']=rev_wip['align_score']
            ret['rev_match']=rev_wip['match_str']
            ret['rev_mis_matches']=rev_wip['mis_matches']
            ret['rev_three_prime_mm']=rev_wip['hit']['three_prime_mm']
            ret['rev_three_prime_overhang']=rev_wip['hit']['three_prime_overhang']
        else:
            ret['rev_align_score']=0
        return ret
    def MMHybTm(self, seqA, seqB, DNAConc = None):
        TmA, TmB, TmAB, deltaTm = self.SWTm_WillItHyb(seqA, seqB, DNAConc = DNAConc)[:4]
        return TmAB
    def SWTm_SequenceSimilarity(self, seqA, seqB, DNAConc = None):
        #Are two sequences in a list similar to each other
        #SeqCompare: SeqA-Sense to SeqB-Sense
        #ThermoCompare: SeqA-Hit[0] to SeqB-Hit[1].Complement
        hit = self.PureSeqCompare(seqA, seqB)
        #ExampleOutput:
        #['GAATTC', 'GAATTC']
        TmA =  self.HybTm(seqA, DNAConc = DNAConc)
        TmB =  self.HybTm(seqB, DNAConc = DNAConc)
        if hit[0]:
            sA = hit[0]
            #GAATTC
            #gaattc - ANTISENSE
            #CTTAAG - TRANSLATE
            sB = hit[1].translate(makeTran)
            TmAB = self.HybTm(sA, seqB = sB, DNAConc = DNAConc)
            deltaTm = abs(min(TmA,TmB)-TmAB)
            return TmA, TmB, TmAB, deltaTm, hit[0], hit[1]
        else:
            return TmA, 0.0, TmA, TmA, '', ''
    def AreTheySimilar(self, seqA, seqB, DNAConc = None, dTm = 15.0):
        TmA, TmB, TmAB, deltaTm = self.SWTm_SequenceSimilarity(seqA, seqB, DNAConc = DNAConc)[:4]
        if deltaTm < dTm:
            return True
        else:
            return False
#    def ThermoSecStructSW(self, seqA, seqB):
#
#        hit = MB.SW(seqA,seqB, match = 1, MM = -1, gapOpen = -2, gapExtend = -1,Verbose = Verbose, DEBUG = 0)
#        #ExampleOutput:
#        #['GATCGTAGCTGATCGTAGCTATCGATGCT', 'GATCGTAGCTGATCGTAGCTATCGATGCT']
#        SS_DNAConc = 1.0 ##secStructure
#        if len(hit[0])>3:
#            sA = hit[0]
#            sB = hit[1].translate(makeTran)
#            
#            TmAB = self.HybTm(sA, seqB = sB, DNAConc = SS_DNAConc)
#            return TmAB, sA, sB
#        else:
#            return 0.0, '', ''
    def MaxThermoSecStruct(self, inSeq , wordSize = 4):
        seqLen = len(inSeq)
        tLst = []
        for i in range(seqLen-wordSize):
            qseq = inSeq[i:i+wordSize]
            testSeqLst = [inSeq[:i], inSeq[i+wordSize:]]
            if Verbose: print(qseq + '  ' + str(testSeqLst))
            for tseq in testSeqLst:
                if len(tseq) >= wordSize:
                    TmA, TmB, TmAB, deltaTm, sA, sB = self.SWTm_WillItHyb(qseq, tseq, DNAConc = 1.0)
                    tLst.append([TmAB, sA, sB])
        tLst.sort(lambda x,y: cmp(x[0], y[0]))
        if Verbose: print(str(tLst))
        if tLst:
            winner = tLst[-1]
        else:
            winner = [0.0, '','']
        return winner
    def SecStructTm(self, inSeq,  windowSize = 12):
        return self.MaxThermoSecStruct(inSeq , wordSize = windowSize)[0]
        




TM = ThermoVals()
        
if __name__ == '__main__':
    Verbose =1

        
    s = 'GATCGTAGCTGATCGTAGCTAT'
    print ('\n\n\nSequenceSimilarity start:' + s)
    r = TM.SWTm_SequenceSimilarity('GATCGTAGCTGATCGTAGCTAT', 'GATCGTAGCTTCGTAGCTAT')
    print('r is :' + str(r) + '\n\n\n')
    t = TM.AreTheySimilar('GATCGTAGCTGATCGTAGCTAT', 'GATCGTAGGGTTCGTAGCTAT', DNAConc = None, dTm = 10.0)
    print('AreTheySimilar: %s' %str(t))
    
    s = 'GATCGTAGCTGATCGTAGCTAT'
    print('\n\n\nWillItHyb start:' + s)
    r = TM.SWTm_WillItHyb('GATCGTAGCTGATCGTAGCTAT', 'ATAGCTACGA-CAGCTACGATC')
    print('r is :' + str(r) + '\n\n\n')
    t = TM.MMHybTm(s, Antisense(s), DNAConc = None)
    print('MMHybTm: %.2f' %t)
    
    seq = 'CGGGGCGGAAGgactgctccgctccgcgttgacatca'
    r = TM.SecStructTm(seq, windowSize = 15)
    print(str(r))
    sTm = MB.HairpinTm(seq, MinAnneal = 4)#This version does not see gaps or mismatches
    print(str(sTm))
