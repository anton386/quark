import sys
from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq

class GenBank(object):


    def __init__(self, genbank=None, email=None):
        
        self.genbank = genbank
        self.email = email
    

    def save(self, gi, email):
        
        Entrez.email = email
        handle = Entrez.efetch(db="nucleotide", id=gi, rettype="gb", retmode="text")

        out = open(self.genbank, "w")
        out.write(handle.read())
        out.close()

        handle.close()


    def load(self):
        self.record = SeqIO.read(self.genbank, "genbank")


    def get_cds_from_genbank(self, first_entry=True):
        
        list_of_cds = []
        
        for feature in self.record.features: # Returns list of SeqFeature
            if feature.type == "CDS":
                list_of_cds.append((feature.location.start, feature.location.end))

        if first_entry:
            return list_of_cds[0]
        else:
            return list_of_cds
        

if __name__ == "__main__":
    gb = GenBank(sys.argv[1])
    gb.load()
    cds = gb.get_cds()
    print gb.record.seq[cds[0]:cds[1]].translate(to_stop=True)

    
