Verbose = None

class TwoD_Array:
    def __init__(self, IFN = None, inData = None):
        #Assumes data type is float
        #NullValue = max
        self.Max = None
        self.Min = None
        self.dict = {}
        if inData:
            self.prepdata(inData)
        elif IFN:
            inData = IFN.inDB()
            self.prepdata(inData)
        else:
            self.indata = []
    def __getitem__(self, val):
        return self.dict[val]
    def getVal(self, firstVal, secondVal):
        firstVal = firstVal.upper()
        secondVal = secondVal.upper()
        if firstVal in self.dict:
            firstDict = self.dict[firstVal]
            if secondVal in firstDict:
                return firstDict[secondVal]
            else:
                raise Exception('%s NOT in second dict' %secondVal)
        else:
            raise Exception('%s NOT in first dict' %firstVal)
    def getVal(self, firstVal, secondVal):
        #Use this to take gaps into account
        firstVal = firstVal.upper()
        secondVal = secondVal.upper()
        Gap = 0
        if '-' in firstVal:
            return self.NullValue/2
        if '-' in secondVal:
            return self.NullValue/2
        if firstVal in self.dict:
            firstDict = self.dict[firstVal]
            if secondVal in firstDict:
                return firstDict[secondVal]
            else:
                return self.NullValue/2
                #raise Exception('%s NOT in second dict' %secondVal)
        else:
            return self.NullValue/2
            raise '%s NOT in first dict' %firstVal
    def setVal(self, firstVal, secondVal, Val):
        firstVal = firstVal.upper()
        secondVal = secondVal.upper()
        if not(firstVal in self.dict):
            self.dict[firstVal] = {}
        self.dict[firstVal][secondVal] = Val
    def prepdata(self, inData):
        self.indata = inData
        self.n = len(self.indata[0]) ##assumes that second keyset is the top line
        secondD_LL = self.indata[0]
        for i in range(1,self.n,1):
            firstD = self.indata[i][0].upper()
            self.dict[firstD] = {}
            for j in range(1, self.n, 1):
                secondD = secondD_LL[j].upper()
                self.dict[firstD][secondD] = self.indata[i][j]
        for i in range(1,self.n,1):
            line = self.indata[i]
            if len(line) != self.n:
                while len(line)<self.n:
                    line.append(None)
            valLst = line[1:]
            for v in range(len(valLst)):
                val = valLst[v]
                if val:
                    valLst[v] = float(val)
                else:
                    valLst[v] = 0
            maxrow = max(valLst)
            minrow = min(valLst)
            self.indata[i] = [line[0]] + valLst
            if self.Max == None:
                self.Max = maxrow
            if self.Min ==None:
                self.Min = minrow
            if maxrow > self.Max:
                self.Max = maxrow
            if minrow > self.Min:
                self.Min = minrow
        self.SetNullValue()
        self.keys = list(self.dict.keys())
        count = 0
        for i in range(len(self.keys)):
            firstkey = self.keys[i]
            for j in range(len(self.keys)):
                secondkey = self.keys[j]
                val = self[firstkey][secondkey]
                if val == '':
                    self[firstkey][secondkey] = self.NullValue
                    count +=1
        if Verbose: print('prepdata: set %i values to default NullValue %.2f' %(count, self.NullValue))
            
    def SetNullValue(self):
        self.NullValue = 2*self.Max
        if Verbose: print('SetNullValue:' + str(self.NullValue))
    def PairDB(self):
        outLst = []
        for i in range(len(self.keys)):
            firstkey = self.keys[i]
            for j in range(len(self.keys)):
                secondkey = self.keys[j]
                outLst.append([firstkey, secondkey, self[firstkey][secondkey]])
        return outLst

        
        
if __name__ == '__main__':
    Verbose = 1
    indata = [['', 'A/T', 'T/A', 'C/G', 'G/C', 'T/G', 'G/T', 'A/G', 'G/A', 'A/C', 'C/A', 'C/T', 'T/C', 'A/A', 'T/T', 'C/C', 'G/G', ],
            ['A/T', '-1.01', '-0.87', '-1.45', '-1.29', '0.07', '0.72', '0.11', '0.01', '0.87', '0.77', '0.64', '0.72', '0.67', '0.65', '1.36', '-0.15', ],
            ['T/A', '-0.59', '-1.01', '-1.31', '-1.46', '0.34', '0.43', '0.48', '0.7', '0.92', '1.33', '0.98', '0.78', '0.7', '0.67', '1.01', '0.48', ],
            ['C/G', '-1.46', '-1.29', '-1.83', '-2.16', '-0.32', '-0.47', '0.01', '0.09', '0.75', '0.79', '0.6', '0.39', '0.4', '-0.1', '0.73', '-0.15', ],
            ['G/C', '-1.31', '-1.45', '-2.23', '-1.83', '-0.59', '0.07', '-0.29', '-0.49', '0.8', '0.48', '0.63', '1.01', '0.14', '0.41', '0.84', '-1.1', ],
            ['T/G', '0.43', '0.72', '0.07', '-0.47', '0.74', '0.52', '', '', '', '', '', '', '', '', '', '', ],
            ['G/T', '0.34', '0.07', '-0.59', '-0.32', '1.15', '0.74', '', '', '', '', '', '', '', '', '', '', ],
            ['A/G', '0.7', '0.01', '-0.49', '0.09', '', '', '', '', '', '', '', '', '', '', '', '', ],
            ['G/A', '0.48', '0.11', '-0.29', '0.01', '', '', '', '', '', '', '', '', '', '', '', '', ],
            ['A/C', '1.33', '0.77', '0.48', '0.79', '', '', '', '', '', '', '', '', '', '', '', '', ],
            ['C/A', '0.92', '0.87', '0.8', '0.75', '', '', '', '', '', '', '', '', '', '', '', '', ],
            ['C/T', '0.78', '0.72', '1.01', '0.39', '', '', '', '', '', '', '', '', '', '', '', '', ],
            ['T/C', '0.98', '0.64', '0.63', '0.6', '', '', '', '', '', '', '', '', '', '', '', '', ],
            ['A/A', '0.7', '0.67', '0.14', '0.4', '', '', '', '', '', '', '', '', '', '', '', '', ],
            ['T/T', '0.67', '0.65', '0.41', '-0.1', '', '', '', '', '', '', '', '', '', '', '', '', ],
            ['C/C', '1.01', '1.36', '0.84', '0.73', '', '', '', '', '', '', '', '', '', '', '', '', ],
            ['G/G', '0.48', '-0.15', '-1.1', '-0.15', '', '', '', '', '', '', '', '', '', '', '', '', ],
                ]
    
    twoD = TwoD_Array(IFN = None, inData = indata)
    pdb = twoD.PairDB()
    for obj in pdb:
        print (str(obj))
    
