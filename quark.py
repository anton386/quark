import sys
import getopt
import collections
import math

import multiprocessing as mp
import functools

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
from variants import Variants
from translation import Translation
from genbank import GenBank

# as long as the end of the window is less than the
#    - total length of the genome
#    - the stated end position
# while (pos + interval) < length_of_genome and (pos + interval) < self.end:
def parallel(pos, sam, con, window_range, increment_window):
                
    win = Window(pos,
                 window_range,
                 increment_window)
    hap = Haplotypes(pysam.Samfile(sam, "rb"), con, win)

    # Calculate the Entropy here
    haplotype_frequency = []
    total_coverage = 0
    total_haplotypes = 0
    for k1, v1 in hap.haplotype.items():
        for k2, v2 in v1.items():
            haplotype_frequency.append(v2)
            total_coverage += v2
            total_haplotypes += 1
        
    haplotype_relative_freq = [float(h)/float(total_coverage) for h in haplotype_frequency]
    entropy = 0.00
    for h in haplotype_relative_freq:
        entropy += h * math.log(h)
    entropy = -1.0 * entropy

    row = (pos, entropy, total_coverage, total_haplotypes)
    sys.stdout.write("\t".join(map(str, list(row))) + "\n")
    return row

class Quark(object):

    def __init__(self, command, filename=None, reference=None, allele_of_interest=None,
                 window_range=None, increment_window=None, use_consensus_as_ref=None,
                 export_to_fasta=None, export_to_vcf=None, format=2, out=None,
                 start=None, end=None, option=None, threads=None, gb=None):

        self.command = command
        self.filename = filename
        self.sam = pysam.Samfile(self.filename, "rb")
        self.ref = Reference(reference)
        self.allele_of_interest = allele_of_interest
        self.window_range = window_range
        self.increment_window = increment_window
        self.use_consensus_as_ref = True if self.command == "consensus" else use_consensus_as_ref
        self.export_to_fasta = export_to_fasta
        self.export_to_vcf = export_to_vcf
        self.format = format
        self.out = out
        self.start = start
        self.end = end
        self.option = option
        self.threads = threads
        self.gb = gb

        # once done, execute
        self._start()    
        self.sam.close()


    def _start(self):

        if self.use_consensus_as_ref:
            con = ConsensusReference(self.filename, self.ref, self.export_to_fasta)
        else:
            con = self.ref

        if self.command == "top_ten":
            win = Window(self.allele_of_interest,
                         self.window_range,
                         self.increment_window)  # object contains all coordinates
            hap = Haplotypes(self.sam, con, win)  # iterate through reads, and then win sizes
            self.top_ten(hap, win)
        elif self.command == "top_one":
            win = Window(self.allele_of_interest,
                         self.window_range,
                         self.increment_window)  # object contains all coordinates
            hap = Haplotypes(self.sam, con, win)  # iterate through reads, and then win sizes
            self.top_one(self.allele_of_interest, hap, win)
        elif self.command == "rank_ten":
            pass
            #self.rank_ten(hap, win, fix_window)
        elif self.command == "haplotype_calling":
            win = Window(self.allele_of_interest,
                         self.window_range,
                         self.increment_window)  # object contains all coordinates
            hap = Haplotypes(self.sam, con, win)  # iterate through reads, and then win sizes
            self.haplotype_calling(hap, win, con, self.format, self.out)  # TODO format
        elif self.command == "windows":
            pass
        elif self.command == "variant_calling":
            region = None
            threshold = None
            max_depth = 200000
            var = Variants(self.filename, con, region, max_depth, threshold)  # iterate through reads, and then win sizes
            self.variant_calling(var)
        elif self.command == "consensus":
            con.export_to_fasta()
        elif self.command == "haplotype_entropy":
            
            # exit if this requirement fails
            if "-" in self.window_range:
                sys.stderr("[Error] Only a single window value is allowed for entropy_calling")
                sys.exit(99)

            # define parameters
            length_of_genome = len(con.sequence)
            interval = int(self.window_range)/2
            pos = self.start + interval

            # Start the multiprocessing module
            pool = mp.Pool(processes=self.threads)
            closure = functools.partial(parallel, sam=self.filename, con=con,
                                        window_range=self.window_range,
                                        increment_window=self.increment_window)
            
            if self.option == "fixed_window":
                end = min(length_of_genome, self.end)
                data = pool.map(closure, range(pos, end, int(self.window_range)))
                pool.close()
                pool.join()
            elif self.option == "moving_window":
                end = min(length_of_genome, self.end)
                data = pool.map(closure, range(self.start, end))
                pool.close()
                pool.join()

            # Send for plotting
            self.haplotype_entropy(data, self.out)

    def variant_calling(self, var):
        #for k1, v1 in sorted(var.variants.items(), key=lambda q: q[0]):
        #    print k1, v1, var.coverage[k1], var.calculate_mean_quality_score(k1)

        var.export_to_vcf()

    def haplotype_calling(self, hap, win, ref, format, out):
        outfile = open(out, "w")

        if format == 0 or format == 1:
            # collect set of variant positions
            unique = []
            for k1, v1 in hap.haplotype.items():
                for k2 in v1.keys():
                    for k3 in k2:
                        unique.append(k3[0])

            sorted_unique = sorted(set(unique))
            if format == 1:
                outfile.write(",".join(["id"] + map(str, sorted_unique) + ["frequency"]) + "\n")

            # collect the bases at this positions
            for k1, v1 in hap.haplotype.items():
                for n, (k2, v2) in enumerate(sorted(v1.items(), key=lambda q: q[1], reverse=True)):
                    out = collections.OrderedDict()

                    for u in sorted_unique:  # set([2816...2819])
                        added = False
                        for k3 in k2:
                            if k3[0] == u:
                                # check whether its a deletion greater than 2
                                if "-" in k3[2]:
                                    if len(k3[2]) > 1:
                                        for i in range(len(k3[2])):
                                            out[u+i] = "-"
                                    else:
                                        out[u] = k3[2]
                                elif "I" in k3[1]:
                                    out[u] = ref.sequence[u-1] + k3[2]
                                else:
                                    out[u] = k3[2]
                                added = True
                        if not added:
                            try:
                                out[u]
                            except KeyError:
                                out[u] = ref.sequence[u-1]

                    if format == 0:
                        outfile.write("\t".join(["".join(out.values()), str(v2)]) + "\n")
                    elif format == 1:
                        outfile.write(",".join([str(n)] + out.values() + [str(v2)]) + "\n")
                        
        elif format == 2:

            for k1, v1 in hap.haplotype.items():
                for n1, (k2, v2) in enumerate(sorted(v1.items(), key=lambda q: q[1], reverse = True)):
                    outfile.write(str(n1) + "\t" + str(k2) + "\t" + str(v2) + "\n")

        elif format == 3:

            # TEST
            if self.gb:
                gb = GenBank(self.gb)
                gb.load()
                cds = gb.get_cds_from_genbank()
            
                t = Translation(ref, cds[0], cds[1])
            
            for k1, v1 in hap.haplotype.items():
                for n1, (k2, v2) in enumerate(sorted(v1.items(), key=lambda q: q[1], reverse = True)):

                    # TEST
                    ref_protein = t.get_translation(ref.sequence, start=win.windows[0][0], end=win.windows[0][1])
                    if self.gb:
                        s = t.load_variants(k2)
                        q = t.get_translation(s, start=win.windows[0][0], end=win.windows[0][1])
                        p = t.diff(ref_protein, q)
                    
                        outfile.write(str(n1) + "\t" + str(k2) + "\t" + str(p) + "\t" + str(v2) + "\n")
                    else:
                        outfile.write(str(n1) + "\t" + str(k2) + "\t" + str(v2) + "\n")
            
        
        outfile.close()

    def haplotype_entropy(self, data, outfile):
        plt.plot([d[0] for d in data], [d[1] for d in data])
        plt.xlabel("Genomic Coordinates")
        plt.ylabel("Entropy")
        plt.title("Entropy across the genome")
        plt.savefig(outfile, dpi=200)
        

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
python quark.py [top_ten|top_one|rank_ten] <options>
  -b --bam <bam>
  -r --ref <reference>
  -w --window <range_of_window>
  -p --position <pos_of_interest>
  -i --increment <increment - default:10>
  -c --use-consensus-as-ref
  -f --export-consensus-to-fasta <filename>
  -v --export-variants-to-vcf <filename>
  -h --help
"""

def usage2():
    print """
python quark.py windows <options>
  -b --bam
  -r --ref
  -w --window
  -p --position
  -i --increment
  -c --use-consensus-as-ref
  -h --help show this help screen
"""

def usage3():
    print """
python quark.py variant_calling <options>
  -b --bam
  -r --ref
  -p --position
  -i --increment
  -c --use-consensus-as-ref
  -v --export-variants-to-vcf
  -h --help show this help screen
"""

def usage4():
    print """
python quark.py consensus <options>
  -b --bam
  -r --ref
  -f --export-consensus-to-fasta
"""

def usage5():
    print """
python quark.py haplotype_calling <options>
  -b --bam <bam>
  -r --ref <reference>
  -w --window <range_of_window>
  -p --position <pos_of_interest>
  -i --increment <increment - default:10>
  -c --use-consensus-as-ref
  -f --export-consensus-to-fasta <filename>
  -v --export-variants-to-vcf <filename>
  -m --format <haplotype=0, csv=1, raw=2>
  -g --genbank <filename>
  -o --out <out-file>
  -h --help
"""

def usage6():
    print """
python quark.py haplotype_entropy <options>
  -b --bam <bam>
  -r --ref <reference>
  -w --window <range_of_window>
  -s --start <start>
  -e --end <end>
  -c --use-consensus-as-ref
  -t --threads <threads>
  -o --out <out-file>
  -z --option <fixed_window|moving_window>
  -h --help
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
    haplotype_calling
    haplotype_entropy
    consensus
    help
"""

if __name__ == "__main__":

    cmd = sys.argv[1]

    # defines which cmd requires which parameter
    if cmd in ("top_ten", "top_one", "rank_ten"):
        try:
            short_options = "b:r:w:p:i:cf:v:h"
            long_options = ["bam", "ref", "window", "increment",
                            "position", "help",
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
            short_options = "b:r:w:p:i:ch"
            long_options = ["bam", "ref", "window", "increment",
                            "position", "use-consensus-as-ref",
                            "help"]
            opts, args = getopt.getopt(sys.argv[2:], short_options, long_options)
        except getopt.GetoptError as err:
            print str(err)
            usage2()
            sys.exit(99)
    elif cmd == "variant_calling":
        try:
            short_options = "b:r:p:i:chv:"
            long_options = ["bam", "ref", "increment", "position",
                            "use-consensus-as-ref",
                            "export-variants-to-vcf",
                            "help"]
            opts, args = getopt.getopt(sys.argv[2:], short_options, long_options)
        except getopt.GetoptError as err:
            print str(err)
            usage3()
            sys.exit(99)
    elif cmd == "consensus":
        try:
            short_options = "b:r:f:h"
            long_options = ["bam", "ref", "export-consensus-to-fasta",
                            "help"]
            opts, args = getopt.getopt(sys.argv[2:], short_options, long_options)
        except getopt.GetoptError as err:
            print str(err)
            usage4()
            sys.exit(99)
    elif cmd == "haplotype_calling":
        try:
            short_options = "b:r:w:p:i:cf:v:o:m:g:h"
            long_options = ["bam", "ref", "window", "increment",
                            "position", "help",
                            "use-consensus-as-ref",
                            "export-consensus-to-fasta",
                            "export-variants-to-vcf",
                            "genbank",
                            "out"]
            opts, args = getopt.getopt(sys.argv[2:], short_options, long_options) 
        except getopt.GetoptError as err:
            print str(err)
            usage5()
            sys.exit(99)
    elif cmd == "haplotype_entropy":
        try:
            short_options = "b:r:s:e:w:co:hz:t:"
            long_options = ["bam", "ref",
                            "start", "end", "window",
                            "use-consensus-as-ref",
                            "out", "help", "option",
                            "threads"]
            opts, args = getopt.getopt(sys.argv[2:], short_options, long_options)
        except getopt.GetoptError as err:
            print str(err)
            usage6()
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
    export_to_fasta = None
    export_to_vcf = None
    format = 2
    out = None
    start = None
    end = None
    option = None
    threads = 1
    gb = None
    
    for o, a in opts:
        if o in ("-b", "--bam"):
            filename = a
        elif o in ("-r", "--ref"):
            reference = a
        elif o in ("-w", "--window"):
            window_range = a
        elif o in ("-i", "--increment"):
            increment = int(a)
        elif o in ("-p", "--position"):
            position = int(a)
        elif o in ("-c", "--use-consensus-as-ref"):
            use_consensus_as_ref = True
        elif o in ("-f", "--export-consensus-to-fasta"):
            export_to_fasta = a
        elif o in ("-v", "--export-variants-to-vcf"):
            export_to_vcf = a
        elif o in ("-m", "--format"):
            format = int(a)
        elif o in ("-o", "--out"):
            out = a
        elif o in ("-s", "--start"):
            start = int(a)
        elif o in ("-e", "--end"):
            end = int(a)
        elif o in ("-z", "--option"):
            option = a
        elif o in ("-t", "--threads"):
            threads = int(a)
        elif o in ("-g", "--genbank"):
            gb = a
        elif o in ("-h", "--help"):
            if cmd in ("top_ten", "top_one", "rank_ten"):
                usage1()
                sys.exit(99)
            elif cmd in ("windows"):
                usage2()
                sys.exit(99)
            elif cmd in ("variant_calling"):
                usage3()
                sys.exit(99)
            elif cmd in ("consensus"):
                usage4()
                sys.exit(99)
            elif cmd in ("haplotype_calling"):
                usage5()
                sys.exit(99)
            elif cmd in ("haplotype_entropy"):
                usage6()
                sys.exit(99)
    

    q = Quark(cmd, filename=filename, reference=reference,
              allele_of_interest=position,
              increment_window=increment,
              window_range=window_range,
              use_consensus_as_ref=use_consensus_as_ref,
              export_to_fasta=export_to_fasta,
              export_to_vcf=export_to_vcf,
              format=format,
              out=out,
              start=start,
              end=end,
              option=option,
              threads=threads,
              gb=gb)
