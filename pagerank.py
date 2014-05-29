import sys
import random
import math
import numpy as np
import networkx as nx
import networkx.algorithms.link_analysis as la
import networkx.algorithms.centrality as centrality

from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq

# Given a Graph data structure,
# perform PageRank to get rank
# the viral haplotypes that are interesting
# Questions:
#   How do I represent the Graph
#   How do I traverse the Graph

class GenBank(object):


    def __init__(self, genbank=None, email=None):
        self.genbank = genbank
        self.email = email


    def save(self, gi):
        Entrez.email = self.email
        handle = Entrez.efetch(db="nucleotide", id=gi, rettype="gb", retmode="text")

        out = open(self.genbank, "w")
        out.write(handle.read())
        out.close()
        
        handle.close()

    
    def load(self):
        self.record = SeqIO.read(self.genbank, "genbank")


    def pprint(self):
        print self.record


class PageRank(object):

    def __init__(self, pg):
        self.pg = pg
        self.simulations = 10000000
        self.p_restart = 0.9

    def rank_it(self):
        rank = la.pagerank(self.pg.graph, max_iter=self.simulations, alpha=self.p_restart)
        
        print "\t".join(["id","pagerank","freq-abs","freq-rel","centrality"])
        for k, v in sorted(rank.items(), key=lambda q: q[1], reverse=True):
            print "%s:\t%s\t%s\t%s\t%s" % (k, v,
                                           self.pg.freq[k],
                                           self.pg.freq[k]/float(self.pg.total),
                                           self.pg.central[k])

class PageGraph(object):
    def __init__(self, csv):
        self.csv = csv
        self.markers = []
        self.samples = 0
        self.freq = {}
        self.total = 0

        self.data = {}
        self.matrix = None
        self.graph = None
        self.central = None
        
    def load_haplotype_data(self):
        with open(self.csv, "r") as f:
            for n, g in enumerate(f):
                h = g.strip("\r\n").split(",")
                if n == 0:
                    self.markers = map(int, h[1:-1])
                else:
                    self.data[int(h[0])] = h[1:-1]
                    self.freq[int(h[0])] = int(h[-1])
                    self.total += int(h[-1])
                    self.samples += 1

    def distance_matrix(self, distance_function="manhanttan"):
        x = []
        for i in range(self.samples):
            y = []
            for j in range(self.samples):
                y.append(0)
            x.append(y)
        
        self.matrix = np.array(x)
        for i in range(self.samples):
            for j in range(i):
                if i != j:
                    distance = 0
                    for var1, var2 in zip(self.data[i], self.data[j]):
                        if var1 != var2:
                            distance += 1

                    self.matrix[i][j] = distance
                    self.matrix[j][i] = distance
                else:
                    self.matrix[i][j] = 0

    def graph_it(self):
        G = nx.Graph()
        G.add_nodes_from(range(len(self.matrix)))

        edges = []
        for i in range(len(self.matrix)):
            for j in range(len(self.matrix[i][:i])):
                if i != j:
                    if self.matrix[i][j] < 2:
                        #edges.append((i, j, self.matrix[i][j]))
                        edges.append((i, j))
        G.add_edges_from(edges)
        self.graph = G

    def central_it(self):
        self.central = centrality.degree_centrality(pg.graph)
        


if __name__ == "__main__":

    #gb = GenBank(genbank=sys.argv[1])
    #gb.load()
    #gb.pprint()
    
    pg = PageGraph(sys.argv[1])
    pg.load_haplotype_data()
    pg.distance_matrix()
    pg.graph_it()
    pg.central_it()

    pr = PageRank(pg)
    pr.rank_it()

    
    
    
