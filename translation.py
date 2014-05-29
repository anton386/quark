import sys
import soma.ref as ref
from Bio.Seq import Seq

class Translation(object):

    def __init__(self, ref, cds_start, cds_end):
        self.ref = ref
        self.cds_start = cds_start
        self.cds_end = cds_end

    def load_variants(self, variants):
        seq = list(self.ref.sequence)
        for pos, typev, variant in variants:
            if "I" in typev:
                seq[pos-1] = variant
            elif "D" in typev:
                for i in range(len(variant)):
                    seq[pos-1+i] = ""
            elif "S" in typev:
                seq[pos-1] = variant
        
        return "".join(seq)

    def load_variants_from_csv(self):
        pass

    def get_translation(self, seq, start=None, end=None):
        s = self.convert_to_aa(start)
        e = self.convert_to_aa(end)
        
        return Seq(seq)[self.cds_start:self.cds_end].translate()[s:e]

    def convert_to_aa(self, pos):
        return (pos - self.cds_start + 1)/3

    def diff(self, seq1, seq2, offset=0):
        variants = []
        for n, (v1, v2) in enumerate(zip(seq1, seq2)):
            if v1 != v2:
                variants.append((n+offset, v2))

        return tuple(variants)
                

if __name__ == "__main__":
    r = ref.Reference(sys.argv[1])
    
    gb = GenBank(sys.argv[2])
    gb.load()
    cds = gb.get_cds_from_genbank()
    
    t = Translation(sys.argv[1], cds[0], cds[1])

    for variant in variants:
        s = t.load_variants()
        print t.get_translation(s, start=0, end=0)
    

'''
if __name__ == "__main__":
    
    r = ref.Reference(sys.argv[2])
    
    if sys.argv[1] == "canon":
        print Seq(r.sequence[94:]).translate(to_stop=True)
    elif sys.argv[1] == "alt1_w_indel":
        # insert indels
        seq1 = list(r.sequence)
        print Seq("".join(seq1[94:2795] + seq1[2800:])).translate(to_stop=True)
        #DEBUG print "".join(seq1[:2795] + ["-" for i in range(5)] + seq1[2800:])
    elif sys.argv[1] == "alt1_w_snp_and_indel":
        # insert snp and indels
        seq1 = list(r.sequence)
        seq1[2790] = "T"
        seq1[2791] = "C"
        seq1[2792] = "C"
        seq1[2793] = "A"
        print Seq("".join(seq1[94:2795] + seq1[2800:])).translate(to_stop=True)
    elif sys.argv[1] == "test":
        seq1 = list(r.sequence)
        if int(sys.argv[3]) == 0:
            print ">0"
        elif int(sys.argv[3]) == 1:
            seq1[9938-1] = "C"
            print ">1"
        elif int(sys.argv[3]) == 2:
            seq1[9933-1] = "C"
            print ">2"
        elif int(sys.argv[3]) == 3:
            seq1[9911-1] = "C"
            print ">3"
        elif int(sys.argv[3]) == 4:
            seq1[9909-1] = "T"
            print ">4"
        elif int(sys.argv[3]) == 5:
            seq1[9950-1] = "C"
            print ">5"
        
        print Seq("".join(seq1[94:])).translate(to_stop=True)
        
    elif sys.argv[1] == "hcv":
        seq1 = list(r.sequence)
        #seq1[8226-1] = "G"
        #seq1[8228-1] = "G"
        seq1[8223-1] = "C"
        
        print Seq("".join(seq1[341:])).translate(to_stop=True)[2629-1:2629+8]
        #print Seq(r.sequence[341:]).translate(to_stop=True)[2629-1:2629+8]
    elif sys.argv[1] == "hcv2":
        seq1 = list(r.sequence)
        #seq1[8823-1] = "C"
        #seq1[8804-1] = "G"
        #seq1[8806-1] = "G"
        #seq1[8859-1] = "C"

        print Seq("".join(seq1[341:])).translate(to_stop=True)[2828-1:2828+8]
        #print Seq(r.sequence[341:]).translate(to_stop=True)[2828-1:2828+8]
'''
