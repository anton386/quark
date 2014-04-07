import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="darkgrid")

from soma.sam import SAM
from soma.ref import Reference

# quark libraries
from window import Window
from consensus import ConsensusReference
from haplotypes import Haplotypes

class Quark(object):

    def __init__(self, filename, reference, allele_of_interest,
                 smallest_window, largest_window,
                 command = "default", increment_window = 10):
        
        self.sam = SAM(filename).content()
        self.ref = Reference(reference)
        self.filename = filename
        self.allele_of_interest = int(allele_of_interest)
        self.smallest_window = int(smallest_window)
        self.largest_window = int(largest_window)
        self.increment_window = int(increment_window)
        self.command = command.strip()
        
        # once done, execute
        self._start()
        

    def _start(self):
        win = Window(self.allele_of_interest,
                     self.smallest_window,
                     self.largest_window,
                     self.increment_window)  # object contains all coordinates
        con = ConsensusReference("/data/chenga/projects/dengue/data/alignment/S15/DenP-S15.sorted.ol.qc.bam", self.ref)
        hap = Haplotypes(self.sam, con, win)  # iterate through reads, and then win sizes

        if self.command == "top_ten":
            self.top_ten(hap, win)
        elif self.command == "top_one":
            self.top_one(self.allele_of_interest, hap, win)
        elif self.command == "rank_ten":
            pass
            #fix_window = 70
            #self.rank_ten(hap, win, fix_window)
        elif self.command == "distribution":
            self.haplotype_distribution(hap, win)


    def haplotype_distribution(self, hap, win):
        for k1, v1 in hap.haplotype.items():
            freq = []
            counter = 0
            total_reads = 0
            sorted_v = sorted(v1.items(), key=lambda q: q[1], reverse=True)
            for k2, v2 in sorted_v:
                total_reads += v2
            for k2, v2 in sorted_v:
                if counter < 20:
                    freq.append(float(v2)/float(total_reads))
                counter += 1

            plt.hist(freq, label=str(k1), alpha=0.4)

        #plt.xlim(-5, 20)
        #plt.ylim(-0.2, 0.2)
        plt.legend()
        plt.savefig("/data/chenga/haplotype_histogram.png", dpi=200)


    def top_one(self, allele_of_interest, hap, win):
        for k1, v1 in hap.haplotype.items():
            counter = 0
            total_reads = 0
            top_1 = []
            sorted_v = sorted(v1.items(), key=lambda q: q[1], reverse=True)
            for k2, v2 in sorted_v:
                total_reads += v2
            for k2, v2 in sorted_v:
                for k3 in k2:
                    if k3[0] == allele_of_interest and counter < 1:
                        top_1.append((k2, v2))
                        counter += 1
            print "\t".join(map(str, [top_1[0][0],
                                      top_1[0][1],
                                      float(top_1[0][1])/float(total_reads)]))

    def top_ten(self, hap, win):
        for k1, v1 in hap.haplotype.items():
            counter = 0
            total_reads = 0
            top_10 = []
            for k2, v2 in sorted(v1.items(), key=lambda q: q[1], reverse=True):
                total_reads += v2
                if counter < 10:
                    top_10.append((k2, v2))
                counter += 1

            output = map(str, [k1,
                               win.get_window_start_by_size(k1),
                               win.get_window_end_by_size(k1),
                               len(v1),
                               total_reads])

            # print output
            print "\t".join(output)
            print "Top 10 Haplotypes"
            print "\n".join([str(k) + "\t" + \
                             str(v) + "\t" + \
                             str(float(v)/float(total_reads)) for k, v in top_10])
            
            
        

if __name__ == "__main__":
    if sys.argv[1] == "-h":
        print "python quark.py <filename> <reference> <pos_of_interest> <smallest_window> <biggest_window> <command>"
    else:
        filename, reference, alle_of_intrst, smllst_win, bggst_win, command = sys.argv[1:]
        q = Quark(filename, reference, alle_of_intrst,
                  smllst_win, bggst_win, command)
