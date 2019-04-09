from Bio import SeqIO
import os
from Bio.Seq import Seq
from Bio import SearchIO
from bin.parseMafft import parseMafft
class checkLTR_bygBLAST():
    def __init__(self, LTRBlastXml,ltr=3):
        self.infileXML = LTRBlastXml
        self.ltr = ltr
        self.main()

    def main(self):
        self.parseBlastXml()

    def alignHSPfragmentToVector(self, HSPfragment):
        query, hit = str(HSPfragment.query.seq), str(HSPfragment.hit.seq)
        start_query = HSPfragment.query_range[0]
        for i, bases in enumerate(query):
            print(i,bases,hit[i])
        return True

    def parseBlastXml(self):
        """
        parse xml after BLAST LTR vs genome
        :return: vector self.outVector containing numbers
        corresponing to the number of hits covering each base
        in query sequences
        """
        count = 0
        file2 = open(self.infileXML)
        s1 = SearchIO.parse(file2, 'blast-xml')
        count_hits = 0
        #  query
        parsed_clusters = []
        for recor in s1:
            cluster = recor.description.split("::")[1]
            if not cluster in parsed_clusters:
                f_name = str(self.ltr) + "LTR_" + recor.id.replace("|","_") + "Cluster{}".format(cluster) + ".fasta"
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
                parsed_clusters.append(cluster)
                if os.stat(f_name).st_size != 0:
                    parseMafft(f_name, ltr=self.ltr, run_mafft=True)
                else:
                    print("ATTENTION! No sequences in file {}".format(f_name))



checkLTR_bygBLAST(
    r'5_Selected_LTRsCluster39240.fasta_vs_sunflower',
ltr=5)