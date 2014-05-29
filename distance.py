import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import hcluster
import networkx as nx
import math


class Distance(object):

    def __init__(self, csv):
        self.csv = csv
        self.data = {}
        self.markers = []
        self.samples = 0
        self.matrix = None
        self.freq = {}
        self.total = 0

        self.load_haplotype_data()
        self.distance_matrix()
        #self.dendrogram()
        #self.visualize()
        self.graph2()

    def load_haplotype_data(self):
        with open(self.csv, "r") as f:
            for n, g in enumerate(f):
                h = g.strip("\r\n").split(",")
                if n == 0:  # extract header
                    self.markers = map(int, h[1:-1])
                else:
                    self.data[int(h[0])] = h[1:-1]
                    self.freq[int(h[0])] = int(h[-1])
                    self.total += int(h[-1])
                    self.samples += 1

    def distance_matrix(self, distance_function="manhattan"):
        self.matrix = np.array([ [0] * self.samples for i in range(self.samples) ])
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

    def dendrogram(self):
        self.linkage = hcluster.linkage(hcluster.squareform(self.matrix), method="complete")

    def visualize(self):
        plt.figure(figsize=(8,160))
        hcluster.dendrogram(self.linkage, orientation = "right",
                            leaf_font_size = 1,
                            leaf_rotation = 15)
        plt.savefig("DenP-S15_illumina_tree.png", dpi=200)

    def graph(self):
        plt.figure(figsize=(20,20))
        
        G = nx.Graph()
        G.add_nodes_from(range(len(self.matrix)))

        edges = []
        for i in range(len(self.matrix)):
            for j in range(len(self.matrix[i][:i])):
                if i != j:
                    if self.matrix[i][j] < 2:
                        edges.append((i, j, self.matrix[i][j]))

        G.add_weighted_edges_from(edges)

        nx.draw(G, with_labels=False)
        plt.savefig(sys.argv[2], dpi=200)

    def graph2(self):
        plt.figure(figsize=(20,20))

        G = nx.Graph()
        G.add_nodes_from(range(len(self.matrix)))

        edges = []
        for i in range(len(self.matrix)):
            for j in range(len(self.matrix[i][:i])):
                if i != j:
                    if self.matrix[i][j] < 2:
                        edges.append((i, j, self.matrix[i][j]))

        G.add_weighted_edges_from(edges)

        pos = nx.spring_layout(G)
        
        freq = self.freq.values()
        freq = [ math.log(float(f))*100.0 if f != 1 else 10.0 for f in freq ]
        
        nx.draw_networkx_edges(G, pos)
        nx.draw_networkx_nodes(G, pos,
                               nodelist=self.freq.keys(),
                               node_size=freq)
        plt.savefig(sys.argv[2], dpi=200)
    
if __name__ == "__main__":
    d = Distance(sys.argv[1])
