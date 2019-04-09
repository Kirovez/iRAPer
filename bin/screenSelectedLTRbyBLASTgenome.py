from Bio import SeqIO
import os
from Bio.Seq import Seq
from Bio import SearchIO

class checkLTR_bygBLAST():
    def __init__(self, LTRBlastXml,ltr=3):
        self.infileXML = LTRBlastXml
        self.ltr = ltr
        self.created_fastas = []
        self.main()

    def main(self):
        self.created_fastas = self.parseBlastXml()

    def alignHSPfragmentToVector(self, HSPfragment):
        query, hit = str(HSPfragment.query.seq), str(HSPfragment.hit.seq)
        start_query = HSPfragment.query_range[0]
        for i, bases in enumerate(query):
            print(i,bases,hit[i])
        return True

    def parseBlastXml(self):
        """
        parse xml after BLAST LTR vs genome,
        :return: created fasta file containing similar LTR sequences
        corresponing to the number of hits covering each base
        in query sequences
        """
        created_fastas = []
        file2 = open(self.infileXML)
        s1 = SearchIO.parse(file2, 'blast-xml')
        for recor in s1:
            f_name = "{0}_LTR::{1}::{2}.fasta".format(self.infileXML, self.ltr, recor.id)
            outFasta = open(f_name, 'w')
            query_length = recor.seq_len
            number_of_hits = 0
            # hit sequence
            for HSP in recor:
                hs = HSP.hsps
                # HSP of the hit sequence
                for u in hs:
                    #print(u.query.seq)
                    if u.evalue <= 1e-5 and abs(u.query_range[1] - u.query_range[0]) * 100 / query_length >= 75 \
                        and u.ident_num * 100 / query_length > 85:
                        u.hit.id += str(number_of_hits)
                        new_s = str(u.hit.seq).replace('-','')
                        u.hit.seq = Seq(new_s)
                        SeqIO.write(u.hit, outFasta, 'fasta')
                        number_of_hits += 1
            outFasta.close()
            print("Number of sequences in file for {}".format(f_name), number_of_hits)
            if os.stat(f_name).st_size != 0:
                created_fastas.append(f_name)
            else:
                print("ATTENTION! No sequences in file {}".format(f_name))
        print("Number of fasta files created:", len(created_fastas))
        return created_fastas



# checkLTR_bygBLAST(
#     r'5_Selected_LTRsCluster39240.fasta_vs_sunflower',
# ltr=5)