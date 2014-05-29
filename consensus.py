import sys
import pysam

class ConsensusReference(object):

    def __init__(self, bam, reference, out="consensus.fa"):
        # FASTA Specifications
        self.id = ""
        self.sequence = ""
        self.bam = bam
        self.ref = reference
        self.out = out

        self.build()
        
    def set_id(self, id):
        self.id = id.lstrip(">")

    def set_sequence(self, seq):
        self.sequence = seq.strip("\r\n\t ").upper()

    def export_to_fasta(self):
        export = open(self.out, "w")
        export.write(">" + self.id + "\n")
        export.write(self.sequence + "\n")

    def build(self):
        samfile = pysam.Samfile(self.bam, "rb")
        consensus_reference = dict()
        for pileup in samfile.pileup(self.ref.id.split()[0]):
            out = [pileup.pos]
            bases = []
            for pileupread in pileup.pileups:
                if pileupread.indel == 0:
                    bases.append(pileupread.alignment.seq[pileupread.qpos])
                elif pileupread.indel > 0:
                    bases.append(str(pileupread.indel) + "I")  # insertion
                elif pileupread.indel < 0:
                    bases.append(str(-1*pileupread.indel) + "D")  # deletion

            total_bases = 0
            sorted_bases = sorted(set(bases))
            for base in sorted_bases:
                counter = bases.count(base)
                total_bases += counter

            freq = 0
            consensus = self.ref.sequence[pileup.pos]  # taking the initial from ref
            for base in sorted(set(bases)):
                counter = bases.count(base)
                if (counter > freq and
                    (base != "N" or
                     base[-1] != "I" or
                     base[-1] != "D")):
                    freq = counter
                    consensus = base

            consensus_reference[pileup.pos] = consensus

        # build the consensus
        self.id = self.ref.id
        for i in range(len(self.ref.sequence)):
            try:
                self.sequence += consensus_reference[i]
            except KeyError:
                self.sequence += self.ref.sequence[i]

                    
