import collections, sys
import Graph

def twin(km):
    #return Seq.reverse_complement(km)
    revcompl = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A'}[B] for B in x][::-1])
    return revcompl(km)

def kmers(seq,k):
    for i in xrange(len(seq)-k+1):
        yield seq[i:i+k]

def fw(km):
    for x in 'ACGT':
        yield km[1:]+x

def bw(km):
    for x in 'ACGT':
        yield x + km[:-1]

"""
def myParse(f,k):
    h = open(f, "rU")
    for i,line in enumerate(h):
        if (i%4 == 1):
            segments = filter(lambda x: len(x)>=k,line.strip().split("N"))
            for s in segments:
                yield s
"""

"""
def myParse(f,k,startAtLine=0,pfn=False,printProgress=False):
    if printProgress:
        print "myParse(f="+str(f)+", k="+str(k)+", startAtLine="+str(startAtLine)+", printProgress="+str(printProgress)+")"
    assert isinstance(startAtLine, int)
    while not startAtLine%4==0:
        startAtLine-=1
        if printProgress:
            print "startAtLine:",startAtLine
    h = open(f, "rU")
    for i,line in enumerate(h,start=startAtLine):   #We may be skipping some lines
        #print i,line
        if (i%4 == 1):
            segments = filter(lambda x: len(x)>=k,line.strip().split("N"))
            for s in segments:
                if printProgress and (i%400==1):
                    print "myParse is returning the segment from line "+str(i)+" from the file "+str(f)
                yield s,i
"""

def build(readsFromFiles,k=31,limit=1):
    d = collections.defaultdict(int)

    for reads in readsFromFiles:
        for seq_s in reads:
            seq_l = seq_s.split('N')
            for seq in seq_l:
                for km in kmers(seq,k):
                    d[km] += 1
                seq = twin(seq)
                for km in kmers(seq,k):
                    d[km] += 1

    d1 = [x for x in d if d[x] <= limit]
    for x in d1:
        del d[x]

    return d

def contig_to_string(c):
    return c[0] + ''.join(x[-1] for x in c[1:])

def get_contig(d,km):
    c_fw = get_contig_forward(d,km)
    
    c_bw = get_contig_forward(d,twin(km))

    if km in fw(c_fw[-1]):
        c = c_fw
    else:
        c = [twin(x) for x in c_bw[-1:0:-1]] + c_fw
    return contig_to_string(c),c
        

def get_contig_forward(d,km):
    c_fw = [km]
    
    while True:
        if sum(x in d for x in fw(c_fw[-1])) != 1:
            break
        
        cand = [x for x in fw(c_fw[-1]) if x in d][0]
        if cand == km or cand == twin(km):
            break # break out of cycles or mobius contigs
        if cand == twin(c_fw[-1]):
            break # break out of hairpins
        
        if sum(x in d for x in bw(cand)) != 1:
            break

        c_fw.append(cand)

    return c_fw

def all_contigs(d,k):
    done = set()
    r = []
    for x in d:
        if x not in done:
            s,c = get_contig(d,x)
            for y in c:
                done.add(y)
                done.add(twin(y))
            r.append(s)
    
    G = {}
    heads = {}
    tails = {}
    for i,x in enumerate(r):
        G[i] = ([],[])
        heads[x[:k]] = (i,'+')
        tails[twin(x[-k:])] = (i,'-')
    
    for i in G:
        x = r[i]
        #y in fw(last)
        for y in fw(x[-k:]):
            if y in heads:
                G[i][0].append(heads[y])
            if y in tails:
                G[i][0].append(tails[y])
        #z in fw(twin(first))
        for z in fw(twin(x[:k])):
            if z in heads:
                G[i][1].append(heads[z])
            if z in tails:
                G[i][1].append(tails[z])

    return G,r

    

def print_GFA(G,cs,k):
    print "H  VN:Z:1.0"
    for i,x in enumerate(cs):
        print "S\t%d\t%s\t*"%(i,x)
        
    for i in G:
        for j,o in G[i][0]:
            print "L\t%d\t+\t%d\t%s\t%dM"%(i,j,o,k-1)
        for j,o in G[i][1]:
            print "L\t%d\t-\t%d\t%s\t%dM"%(i,j,o,k-1)
    
def createGraphObject(G,cs,k,GraphObject,pfn=False,ps=False):
    if pfn:
        print createGraphObject.__name__
    #Add all the contigs to the graph
    for cID,c in enumerate(cs):
        GraphObject.contigs[cID] = [c, [],[],0]
    #GraphObject.printContigs()

    #Set the IN/OUT of each contig
    #cID is the ID of contig c
    if ps:
        print "About to iterate through G"
    for c_ID in G:
        if ps:
            print c_ID, GraphObject.contigs[c_ID][0],G[c_ID]
            print "inside the for loop"
            print "c_ID:", c_ID
            print "G[c_ID]:", G[c_ID]

        #def connect_a_to_b_ifNotAlreadyConnected(self,aID,bID,A,B):
        #iterate through the IDs of every contig N that c connects to
        for N_ID,plusMinus in G[c_ID][0]:
            if ps:
                print "inside the first for loop"
            if plusMinus == "+":
                #c->N
                GraphObject.connect_a_to_b_ifNotAlreadyConnected(c_ID,N_ID,True,True)
                #GraphObject.addOUT(c_ID,(N_ID,True))
                #GraphObject.addIN(N_ID,(c_ID,True))
            elif plusMinus == "-":
                #c->twin(N)
                #N->twin(c)
                GraphObject.connect_a_to_b_ifNotAlreadyConnected(c_ID,N_ID,True,False)
                #GraphObject.addOUT(c_ID,(N_ID,False))
                #GraphObject.addOUT(N_ID,(c_ID,False))
            else:
                raise Exception('illegal input')

        #iterate through the IDs of every contig N that twin(c) connects to
        for N_ID,plusMinus in G[c_ID][1]:
            if ps:
                print "inside the second for loop"
                print N_ID,plusMinus
            if plusMinus == "+":
                #twin(c)->N
                #twin(N)->c
                GraphObject.connect_a_to_b_ifNotAlreadyConnected(c_ID,N_ID,False,True)
                #GraphObject.addIN(c_ID,(N_ID,False))
                #GraphObject.addIN(N_ID,(c_ID,False))
            elif plusMinus == "-":
                if c_ID==N_ID:
                    continue
                #twin(c)->twin(N)
                #N->c
                GraphObject.connect_a_to_b_ifNotAlreadyConnected(c_ID,N_ID,False,False)
                #GraphObject.addIN(c_ID,(N_ID,True))
                #GraphObject.addOUT(N_ID,(c_ID,True))
            else:
                raise Exception('illegal input')
        
    #Now we need to fix the inner variables of the GraphObject.
    #I.e. update self.ID and self.kmers

    #Find the highest ID used:
    maxID = -1
    for ID in GraphObject.contigs:
        if ID>maxID:
            maxID = ID

    #update the ID:
    if not maxID==-1:   #i.e. empty Graph
        i = GraphObject.getID()
        while i < maxID:
            i = GraphObject.getID()

    #Add the kmers to the dict
    for ID in GraphObject.contigs:
        GraphObject.addKmersFromContig(ID)


if __name__ == "__main__":
    k = int(sys.argv[1])
    fileReads = []
    for f in sys.argv[2:]:
        fileReads.append(myParse(f))    #index i stores a list of all reads from file i
    d = build(fileReads,k,1)
    G,cs = all_contigs(d,k)
    print_GFA(G,cs,k)
    GraphObject = createGraphObject(G,cs,k)
    GraphObject.printContigs()    