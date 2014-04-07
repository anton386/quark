import re
import collections

class Haplotypes(object):

    def __init__(self, sam_obj, ref_obj, window_obj):
        self.sam = sam_obj
        self.ref = ref_obj
        self.win = window_obj
        self.regex_cigar_num_base_to_skip = re.compile("[0-9]+")
        self.regex_cigar_alpha_operator = re.compile("[A-Z]")
        self.get_haplotype_for_each_window()

    def get_haplotype_for_each_window(self):
        self.haplotype = collections.OrderedDict()
        
        for samy in self.sam:  # iterate through each read
            # is our allele found in this read?
            start, end = (int(samy.pos), int(samy.pos) + self.calculate_read_length(samy) - 1)
            if (self.win.allele_of_interest >= start and
                self.win.allele_of_interest <= end):
                # for every window, check whether the span matches the read
                for win in self.win.windows:
                    # check window
                    if (win[0] >= start and
                        win[1] <= end):                        
                        size = win[1] - win[0]
                        variants = self.get_variants_in_haplotype(samy, win)
                        try:
                            self.haplotype[size]
                        except KeyError:
                            self.haplotype[size] = {}

                        try:
                            self.haplotype[size][tuple(variants)] += 1
                        except KeyError:
                            self.haplotype[size][tuple(variants)] = 1

    def calculate_read_length(self, samy):
        numeric = map(int, self.regex_cigar_num_base_to_skip.findall(samy.cigar))
        operator = self.regex_cigar_alpha_operator.findall(samy.cigar)

        length = 0
        for c, op in zip(numeric, operator):
            if op == "M":
                length += c
        return length

    def get_variants_in_haplotype(self, samy, win):
        # samy - a SAM read
        # win  - tuple of (start, end)

        # ignore soft clips
        # let's get the insertion from cigar
        # then, let's get the deletion from cigar
        # now, let's get the substitutions
        variants = []
        numeric = map(int, self.regex_cigar_num_base_to_skip.findall(samy.cigar))
        operator = self.regex_cigar_alpha_operator.findall(samy.cigar)

        start = 0
        ref_start = 0
        for c, op in zip(numeric, operator):
            if op == "I":
                position = int(samy.pos) + ref_start
                if (position >= win[0] and
                    position <= win[1]):
                    variants.append((position,
                                     str(c) + str(op),
                                     samy.seq[start:start+c]))
                start += c
            elif op == "D":
                position = int(samy.pos) + ref_start
                if (position >= win[0] and
                    position <= win[1]):
                    variants.append((position,
                                     str(c) + str(op),
                                     c * "-"))
                ref_start += c
            elif op == "S":
                start += c
            else:  # deal with substitutions
                #print c
                #print start
                #print ref_start
                #print int(samy.pos)+ref_start-1
                #print int(samy.pos)+ref_start-1+c
                #print self.ref.sequence[int(samy.pos)+ref_start-1:int(samy.pos)+ref_start-1+c]
                #print samy.seq[start:start+c]
                for n, (r, s) in enumerate(zip(self.ref.sequence[int(samy.pos)+ref_start-1:int(samy.pos)+ref_start-1+c],
                                               samy.seq[start:start+c])):
                    if r != s:
                        position = int(samy.pos) + ref_start + n
                        if (position >= win[0] and
                            position <= win[1]):
                            variants.append((position,
                                             "S", s))
                ref_start += c
                start += c

        #print samy.qname
        #print variants
        return variants
                        
        
