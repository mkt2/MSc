import Graph

if __name__ == '__main__':
    G = Graph.Graph(k=31,pfn=False,ps=False,al=False,pil=False)
    G.createGraphFromFile("Output/myOutFolder/G_from_BF_counter.txt")
    G.printContigs("The G from file")