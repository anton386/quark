import re
import sys
import pysam
import collections

class Variants(object):

    def __init__(self, bam, ref_obj, region, max_depth=200000, threshold=None):
        self.sam = pysam.Samfile(bam, "rb")
        self.ref = ref_obj
        self.threshold = threshold
        self.max_depth = max_depth
        self.regex_cigar_num_base_to_skip = re.compile("[0-9]+")
        self.regex_cigar_alpha_operator = re.compile("[A-Z]")

        if region:
            self.region_start = int(region.split("-")[0])
            self.region_end = int(region.split("-")[1])
        else:
            self.region_start = None
            self.region_end = None
        
        self.get_variants_from_pileup()


    # fetch every single read version
    def get_variants(self):
        self.variants = collections.OrderedDict()

        for samy in self.sam.fetch(self.ref.id,
                                   start=self.region_start,
                                   end=self.region_end):
            
            variants = self.get_variants_from_read(samy)
            
            for v1 in variants:
                try:
                    self.variants[v1[0]]
                except KeyError:
                    self.variants[v1[0]] = {}

                try:
                    self.variants[v1[0]][(v1[1], v1[2])] += 1
                except KeyError:
                    self.variants[v1[0]][(v1[1], v1[2])] = 1


    def get_mean_quality_and_cov_from_pileup(self):
        pass


    def get_variants_from_read(self, samy):
        variants = []
        numeric = map(int, self.regex_cigar_num_base_to_skip.findall(samy.cigarstring))
        operator = self.regex_cigar_alpha_operator.findall(samy.cigarstring)

        start = 0
        ref_start = 0
        for c, op in zip(numeric, operator):
            if op == "I":
                position = int(samy.pos+1) + ref_start
                # don't add insertions with N
                if "N" not in samy.seq[start:start+c]:
                    variants.append((position,
                                     str(c) + str(op),
                                     samy.seq[start:start+c]))
                    
                start += c
            elif op == "D":
                position = int(samy.pos+1) + ref_start
                variants.append((position,
                                 str(c) + str(op),
                                 c * "-"))
                ref_start += c
            elif op == "S":
                start += c
            else:
                for n, (r, s) in enumerate(zip(self.ref.sequence[int(samy.pos+1)+ref_start-1:int(samy.pos+1)+ref_start-1+c],
                                               samy.seq[start:start+c])):
                    if r != s and s != "N":
                        position = int(samy.pos+1) + ref_start + n
                        variants.append((position,
                                         "S", s))
                ref_start += c
                start += c
        
        return variants
                
            
    # use pileup version
    def get_variants_from_pileup(self):
        self.variants = collections.OrderedDict()
        self.coverage = collections.OrderedDict()
        self.qualities = collections.OrderedDict()
        
        for pileup in self.sam.pileup(self.ref.id, max_depth=self.max_depth):
            
            ref_base = self.ref.sequence[pileup.pos]
            bases = []
            quals = []

            total_bases = 0
            for pileupread in pileup.pileups:
                
                if pileupread.indel == 0:
                    base = pileupread.alignment.seq[pileupread.qpos]
                    if base != "N":
                        if base != ref_base:
                            bases.append(pileupread.alignment.seq[pileupread.qpos])
                            quals.append(ord(pileupread.alignment.qual[pileupread.qpos]) - 33)
                        else:
                            total_bases += 1 # Add all counts of ref bases first
                            quals.append(ord(pileupread.alignment.qual[pileupread.qpos]) - 33)
                        
                elif pileupread.indel > 0:
                    base = pileupread.alignment.seq[pileupread.qpos:pileupread.qpos+pileupread.indel]
                    if "N" not in base:
                        bases.append(base)
                        quals.append(45)
                
                elif pileupread.indel < 0:
                    bases.append(-1*pileupread.indel * "-")
                    quals.append(45)

            sorted_bases = sorted(set(bases))  # get unique base types
            
            for base in sorted_bases:
                counter = bases.count(base)  # count freq. of unique base types
                total_bases += counter  # count total bases

            for base in sorted_bases:
                counter = bases.count(base)
                a_freq = float(counter) / float(total_bases)            
                try:
                    self.variants[pileup.pos + 1][base] = (counter, a_freq)
                except KeyError:
                    self.variants[pileup.pos + 1] = {base: (counter, a_freq)}
            
            for qual in quals:
                try:
                    self.qualities[pileup.pos + 1]
                except KeyError:
                    self.qualities[pileup.pos + 1] = {}
                
                try:
                    self.qualities[pileup.pos + 1][qual] += 1
                except KeyError:
                    self.qualities[pileup.pos + 1][qual] = 1
            
            self.coverage[pileup.pos + 1] = total_bases

    def calculate_mean_quality_score(self, pos):
        counter = 0
        total = 0
        for k, v in self.qualities[pos].items():
            counter += v
            total += k * v
        return float(total)/float(counter)

    def export_to_vcf(self):
        # print header
        print "##fileformat=VCFv4.1"
        print "#%s" % ("\t".join(["CHROM", "POS", "ID",
                                  "REF", "ALT", "QUAL",
                                  "FILTER", "INFO"]), )

        # print content
        for k1, v1 in sorted(self.variants.items(), key=lambda q: q[0]):
            print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (self.ref.id, k1, ".",
                                                      self.ref.sequence[k1],
                                                      self.ret_vcf_variants(v1),
                                                      self.ret_vcf_qual(),
                                                      self.ret_vcf_filter(),
                                                      self.ret_vcf_info(k1, v1))
    def ret_vcf_variants(self, dvar):
        return ",".join(map(str, sorted(dvar.keys(), key=lambda q: q[0])))
    
    def ret_vcf_qual(self):
        return "."

    def ret_vcf_filter(self):
        return "."

    def ret_vcf_info(self, pos, dvar):
        # deal with af
        # deal with ac
        af = []
        ac = []
        for k1, v1 in sorted(dvar.items(), key=lambda q: q[0]):
            ac.append(v1[0])
            af.append(self.convert_exponent(v1[1]))
        af = ",".join(map(str, af))
        ac = ",".join(map(str, ac))
        return "DP=%s;AF=%s;AC=%s;" % (self.coverage[pos],
                                       af, ac)

    def convert_exponent(self, freq):
        return "%.3e" % (freq, )

