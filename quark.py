import sys
import getopt

# other 3rd party libraries
import pysam
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

    def __init__(self, command, filename=None, reference=None, allele_of_interest=None,
                 window_range=None, increment_window=None, use_consensus_as_ref=None,
                 export_to_fasta=None, export_to_vcf=None):

        try:
            self.command = command
            self.filename = filename
            self.sam = pysam.Samfile(self.filename, "rb")
            self.ref = Reference(reference)
            self.allele_of_interest = int(allele_of_interest)
            self.window_range = window_range
            self.increment_window = int(increment_window)
            self.use_consensus_as_ref = True
            self.export_to_fasta = export_to_fasta
            self.export_to_vcf = export_to_vcf

            # once done, execute
            self._start()
            
        except:
            self.sam.close()


    def _start(self):
        win = Window(self.allele_of_interest,
                     self.window_range,
                     self.increment_window)  # object contains all coordinates
        con = ConsensusReference(self.filename, self.ref)
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
            
            
def usage1():
    print """
python quark.py [top_ten|top_one|rank_ten|windows|variant_calling|haplotype_distribution] <options>
  -b --bam <bam>
  -r --ref <reference>
  -w --window <range_of_window>
  -p --position <pos_of_interest>
  -i --increment <increment - default:10>
  -c --use-consensus-as-ref
  -f --export-consensus-to-fasta
  -v --export-variants-to-vcf
  -h --help
"""

def usage2():
    print """
python quary.py windows <options>
  -b --bam
  -r --ref
  -w --window
  -p --interest
  -i --increment
  -c --use-consensus-as-ref
  -h --help show this help screen
"""

def usage3():
    print """
python quary.py variant_calling <options>
  -b --bam
  -r --ref
  -p --interest
  -i --increment
  -c --use-consensus-as-ref
  -v --export-variants-to-vcf
  -h --help show this help screen
"""

def usage_commands():
    print """
python quark.py <cmd> <options>
  list of commands:
    top_ten
    top_one
    rank_ten
    windows
    variant_calling
    haplotype_distribution
    help
"""

if __name__ == "__main__":

    cmd = sys.argv[1]

    # defines which cmd requires which parameter
    if cmd in ("top_ten", "top_one", "rank_ten", "haplotype_distribution"):
        try:
            short_options = "b:r:w:p:icfvh"
            long_options = ["bam", "ref", "window", "increment",
                            "interest", "help",
                            "use-consensus-as-ref",
                            "export-consensus-to-fasta",
                            "export-variants-to-vcf"]
            opts, args = getopt.getopt(sys.argv[2:], short_options, long_options) 
        except getopt.GetoptError as err:
            print str(err)
            usage1()
            sys.exit(99)
    elif cmd == "windows":
        try:
            short_options = "b:r:w:p:ich"
            long_options = ["bam", "ref", "window", "increment",
                            "interest", "use-consensus-as-ref",
                            "help"]
            opts, args = getopt.getopt(sys.argv[2:], short_options, long_options)
        except getopt.GetoptError as err:
            print str(err)
            usage2()
            sys.exit(99)
    elif cmd == "variant_calling":
        try:
            short_options = "b:r:p:ichv"
            long_options = ["bam", "ref", "increment", "interest",
                            "use-consensus-as-ref",
                            "export-variants-to-vcf",
                            "help"]
            opts, args = getopt.getopt(sys.argv[2:], short_options, long_options)
        except getopt.GetoptError as err:
            print str(err)
            usage3()
            sys.exit(99)
    elif (cmd == "help"
          or cmd == "-h"
          or cmd == "--help"):
        usage_commands()
        sys.exit(99)

    # parse specifications to variables here
    filename = None
    reference = None
    window_range = None
    increment = 10
    position = None
    use_consensus_as_ref = False
    export_to_fasta = False
    export_to_vcf = False
    for o, a in opts:
        if o in ("-b", "--bam"):
            filename = a
        elif o in ("-r", "--ref"):
            reference = r
        elif o in ("-w", "--window"):
            window_range = a
        elif o in ("-i", "--increment"):
            increment = int(a)
        elif o in ("-p", "--position"):
            position = int(a)
        elif o in ("-c", "--use-consensus-as-ref"):
            use_consensus_as_ref = True
        elif o in ("-f", "--export-consensus-to-fasta"):
            export_to_fasta = True
        elif o in ("-v", "--export-variants-to-vcf"):
            export_to_vcf = True
        elif o in ("-h", "--help"):
            if cmd in ("top_ten", "top_one", "rank_ten", "distribution"):
                usage1()
                sys.exit(99)
            elif cmd in ("windows"):
                usage2()
                sys.exit(99)
            elif cmd in ("variant_calling"):
                usage3()
                sys.exit(99)

    q = Quark(cmd, filename=filename, reference=reference,
              allele_of_interest=position,
              increment_window=increment,
              window_range=window_range,
              use_consensus_as_ref=use_consensus_as_ref,
              export_to_fasta=export_to_fasta,
              export_to_vcf=export_to_vcf)
