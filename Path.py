import helpers

class Path:
    #contigs is a list of contig tuples
    #L is an integer storing the sum of the cLen of all the contigs in this path
    #   contigs = [(c_ID,True),(x_ID,x_B)]
    #   L = cLen(c,k)+cLen(x,k)
    def __init__(self):
        self.contigs = []
        self.cLengths = []
        self.L = 0

    def __len__(self):
		return self.L

    def lengthOfAllButLast(self):
        assert(len(self.contigs)>1), "The path must contain at least two contigs"
        sumLength = 0
        for i in range(0, len(self.cLengths)-1):
            sumLength += self.cLengths[i]
        return sumLength

    def lenOfAllButFirst(self):
        assert(len(self.contigs)>1), "The path must contain at least two contigs"
        sumLength = 0
        for i in range(1, len(self.cLengths)):
            sumLength += self.cLengths[i]
        return sumLength

    def __contains__(self,c_ID):
        return self.contains(c_ID)

    def isEmpty(self):
        return len(self.contigs)==0

    def getFirstContigTuple(self):
        return self.contigs[0]

    def getLastContigTuple(self):
        return self.contigs[-1]

    def numContigs(self):
        return len(self.contigs)

    def contains(self,c_ID):
        for x_ID,_ in self.contigs:
            if x_ID == c_ID:
                return True
        return False

    def append(self,c_ID,c_B,c,k):
        assert(not c_ID in self), "We can't add a contig twice to the same path"
        self.contigs.append( (c_ID,c_B) )
        l = helpers.cLen(c,k)
        self.cLengths.append(l)
        self.L += l

    def prepend(self,c_ID,c_B,c,k):
        assert(len(self.contigs)>0), "We never use this on an empty path"
        assert(not c_ID in self), "We can't add a contig twice to the same path"
        self.contigs.insert( 0 , (c_ID,c_B) )
        l = helpers.cLen(c,k)
        self.cLengths.insert( 0, l )
        self.L += l

    def equals(self,otherPath):
        if self.L!=otherPath.L:
            return False
        return set(self.contigs)==set(otherPath.contigs)

    def __str__(self):
        return "Path with length: "+str(self.L)+". Contigs: "+str(self.contigs)+". cLengths: "+str(self.cLengths)

    def getContigIDs(self):
        contigIDs = set()
        for x_ID,x_B in self.contigs:
            contigIDs.add(x_ID)
        return contigIDs
    
    def allButFirst(self):
        assert(len(self)>=2), "the path must contain at least 2 contigs"
        newPath = Path()
        newPath.contigs = self.contigs[1:]
        newPath.cLengths = self.cLengths[1:]
        newPath.L = self.L - self.cLengths[0]
        return newPath
    
    def allButLast(self):
        assert(len(self)>=2), "the path must contain at least 2 contigs"
        newPath = Path()
        newPath.contigs = self.contigs[0:-1]
        newPath.cLengths = self.cLengths[0:-1]
        newPath.L = self.L - self.cLengths[-1]
        return newPath

    #Changes the Path into [front_ID, intNodes1, end_ID]
    #Returns the list [front_ID, intNodes1, end_ID]
    #where front_ID and end_ID are integers and intNodes1 is a Path object
    #If front_ID=-1 we take it from the front of the path
    def toFrontIntEnd(self,front_ID=-1):
        if front_ID==-1:
            assert(len(self.contigs)>=3), "otherwise we can't split the path into front_ID, intNodes1, end_ID"
            front_ID = self.contigs[0]
            intNodes1 = Path()
            intNodes1.contigs = self.contigs[1:-1]
            intNodes1.cLengths = self.cLengths[1:-1]
            intNodes1.L = self.L - self.cLengths[0] - self.cLengths[-1]
            end_ID = self.contigs[-1][0]
            return [front_ID, intNodes1, end_ID]
        else:
            assert(len(self.contigs)>=2), "otherwise we can't split the path into intNodes1, end_ID"
            intNodes1 = Path()
            intNodes1.contigs = self.contigs[0:-1]
            intNodes1.cLengths = self.cLengths[0:-1]
            intNodes1.L = self.L - self.cLengths[-1]
            end_ID = self.contigs[-1][0]
            return [front_ID, intNodes1, end_ID]

    #Creates all subsets of this path starting at the first contig
    #but still containing 3 or more contigs (including this path)
    def all3(self):
        l = len(self.contigs)
        assert(l>=3), "otherwise we can't split the path into front_ID, intNodes1, end_ID"
        if l==3:
            return [self]
        return [self]+self.allButLast().all3()

    def all2(self):
        l = len(self.contigs)
        assert(l>=2), "otherwise we can't split the path into intNodes1, end_ID"
        if l==2:
            return [self]
        return [self]+self.allButLast().all2()

class PathList:
    def __init__(self,initPaths=[]):
        self.paths = []
        self.longestPath = 0
        for p in initPaths:
            self.append(p)

    def __len__(self):
		return len(self.paths)

    def append(self,p):
        self.paths.append(p)
        self.longestPath = max( self.longestPath , len(p) )

    def __add__(self, other):
        newPathList = PathList()
        newPathList.paths = self.paths+other.paths
        newPathList.longestPath = max(self.longestPath , other.longestPath)
        return newPathList

    #def appendAllPathsFrom(self, otherPathList):
    #    for p in otherPathList.paths:
    #        self.append(p)

    def __contains__(self,p):
        for x in self.paths:
            if x.equals(p):
                return True
        return False

    def isEmpty(self):
        return len(self.paths)==0

    def pop(self):
        p = self.paths.pop()
        return p

    def getPaths(self):
        return self.paths

    def getContigIDs(self):
        contigIDs = set()
        for p in self.paths:
            for ID in p.getContigIDs():
                contigIDs.add(ID)
        return contigIDs

    def __str__(self):
        if self.isEmpty():
            s = "Empty PathList"
        else:
            s = "PathList:"
            for p in self.paths:
                s += "\n"+str(p)
        return s

    #Returns the length of the longest path currently in the PathList
    def getLongestPath(self):
        return self.longestPath

    def equals(self,other):
        if len(self)!=len(other):
            return False
        for p in self.paths:
            if not p in other:
                return False
        return True

    #Returns a list containing each path in self.paths as [front_ID,intNodes1,end_ID], i.e. a list of lists
    def toListOf_frontIntEnd(self, front_ID=-1):
        potentialBubbles = []
        for p in self.paths:
            front_ID_intNodes1_end_ID = p.toFrontIntEnd(front_ID)
            potentialBubbles.append(front_ID_intNodes1_end_ID)
        return potentialBubbles

    #Returns a new PathList containing p.all3() instead of each path p
    def all3(self):
        newPathList = PathList()
        for p in self.paths:
            p3 = p.all3()
            for p_3 in p3:
                newPathList.append(p_3)
        return newPathList

    def all2(self):
        newPathList = PathList()
        for p in self.paths:
            p2 = p.all2()
            for p_2 in p2:
                newPathList.append(p_2)
        return newPathList

    def filterOutPathsOfLengthLessThan2(self):
        out = PathList()
        for p in self.paths:
            if p.numContigs()>=2:
                out.append(p)
        return out


if __name__ == "__main__":
    p1 = Path()
    p1.append(0,True,"AAAAA",5)
    print p1
    p2 = Path()
    p2.append(0,True,"AAAAA",5)
    P1 = PathList()
    P1.append(p1)
    P2 = PathList()
    P2.append(p2)
    print p1.equals(p2)
    print P1.equals(P2)
    P1.printPaths()
    P2.printPaths()
    print p1 in P1
    print p1 in P2
    print p2 in P1
    print p2 in P2
    print P1
    print p1
    P = PathList([p1])
    print P
    p3 = Path()
    p3.append(0,True,"AAAAA",5)
    p3.append(1,True,"CCCCC",5)
    p3.append(2,True,"TTTTT",5)
    p3.append(3,True,"GGGGG",5)
    p3.append(4,True,"AAAAC",5)
    print "here we go"
    print p3
    allP = p3.all3()
    for p in allP:
        print p
    print "here we go again"
    P3 = PathList([p3])
    print P3
    allP = P3.all3()
    print allP


