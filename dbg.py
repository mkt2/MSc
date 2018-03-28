import collections, sys
import Graph

alphabet = ["A","C","G","T","N"]

def twin(km):
    #return Seq.reverse_complement(km)
    #print km
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

def build(fn,k=31,limit=1):
    d = collections.defaultdict(int)

    for f in fn:
        h = open(f, "rU")
        for lineNr,line in enumerate(h):
            if (lineNr%4 == 1):
                segments = filter(lambda x: len(x)>=k and all([c in alphabet for c in x]),line.strip().split("N"))
                for seq in segments:
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
    print "H\tVN:Z:1.0"
    for i,x in enumerate(cs):
        print "S\t%d\t%s\t*"%(i,x)
        
    for i in G:
        for j,o in G[i][0]:
            print "L\t%d\t+\t%d\t%s\t%dM"%(i,j,o,k-1)
        for j,o in G[i][1]:
            print "L\t%d\t-\t%d\t%s\t%dM"%(i,j,o,k-1)
    
def print_GFA_to_file(G,cs,k,fileName):
    print "print_GFA_to_file(fileName="+str(fileName)+")"
    f = open(fileName, 'w')
    f.write("H\tVN:Z:1.0\n")
    for i,x in enumerate(cs):
        f.write("S\t%d\t%s\t*\n"%(i,x))
        
    for i in G:
        for j,o in G[i][0]:
            f.write("L\t%d\t+\t%d\t%s\t%dM\n"%(i,j,o,k-1))
        for j,o in G[i][1]:
            f.write("L\t%d\t-\t%d\t%s\t%dM\n"%(i,j,o,k-1))
    f.close()

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
    d = build(sys.argv[2:],k,1)
    G,cs = all_contigs(d,k)
    print_GFA(G,cs,k)
    GraphObject = Graph.Graph(k,al=False)
    createGraphObject(G,cs,k,GraphObject)
    GraphObject.saveAs_GFA_toFile("testGFA.gfa")
    #GraphObject.printContigs()    