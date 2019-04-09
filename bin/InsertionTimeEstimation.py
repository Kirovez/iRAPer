from Bio import SeqIO
import os
from Bio.Align.Applications import ClustalwCommandline
class LTR_InsertionTimeCalculator():
    """
    This class takes two fasta files (5'-utr and 3'-utr from LTRdigest output) and perform alignments following by isertion time caluclation based on the formula
    K/2r where:
        K - is distance between two LTRs estimated by Kimura 2-Parameter method K = -0.5 log( (1 - 2p -q) * sqrt( 1 - 2q ) )
        where:
            p = transition frequency
            q = transversion frequency
    r - mutation rate 1.3*10-8 (Ma and Bennetzen, 2004)
    """

    def __init__(self, fasta_LTR_left, fasta_LTR_right, outFile_tab):
        self.fasta_LTR_left = fasta_LTR_left
        self.fasta_LTR_right = SeqIO.index(fasta_LTR_right, "fasta")
        self.r_parameter = 9.4E-09
        self.clustalW2 = "clustalw"
        self.outFile_tab = outFile_tab
        self.run()

    def run(self):
        cnt_l = 0
        cnd_paired = 0
        print("r", self.r_parameter)
        with open(self.outFile_tab, "w") as outfile:
            for seq1 in SeqIO.parse(self.fasta_LTR_left, "fasta"):
                cnt_l += 1
                #print("Sequence ", cnt_l)
                if seq1.id in self.fasta_LTR_right:
                    cnd_paired +=1
                    seq2 = self.fasta_LTR_right[seq1.id]
                    insertion_time = self.align2sequnces(seq1, seq2) # it will align and parse the aligned sequences to insertion time calculation
                    outfile.write(seq2.id + "\t" + str(insertion_time) + "\n")

            print("Number of sequences in left file: ", cnt_l)
            print("Number of sequences in right file: ", len(self.fasta_LTR_right))
            print("Number of valid pairs found: ", cnd_paired)

    def align2sequnces(self,seq1, seq2):
        insertion_time = 0
        # write temporary fasta files with two LTR sequences
        with open("tmp_file.fasta", "w") as fasta_tmp:
            seq1.id += "*1eft"
            SeqIO.write(seq1,fasta_tmp,"fasta")
            SeqIO.write(seq2, fasta_tmp, "fasta")
        #run clustalw
        clustalw_exe = self.clustalW2
        clustalw_cline = ClustalwCommandline(clustalw_exe, infile="tmp_file.fasta")
        # assert os.path.isfile(clustalw_exe), "Clustal W executable missing"
        stdout, stderr = clustalw_cline()

        # open alignment file and count time of insertion
        with open("tmp_file.aln", "rU") as handle:
            seqs =[]
            for record in SeqIO.parse(handle, "clustal"):
                seqs.append(str(record.seq))
                # print(record.seq)
            insertion_time = (self.K2Pdistance(seqs[0], seqs[1]))/(self.r_parameter * 2)
        return (int(insertion_time))

    def K2Pdistance(self, seq1, seq2):
        seq1, seq2 = seq1.upper(), seq2.upper()
        """
        function was borrowed from https://github.com/kgori/python_tools_on_github/blob/master/pairwise_distances.py
        Kimura 2-Parameter distance = -0.5 log( (1 - 2p -q) * sqrt( 1 - 2q ) )
        where:
        p = transition frequency
        q = transversion frequency
        """
        from math import log, sqrt
        pairs = []

        # collect ungapped pairs
        for x in zip(seq1, seq2):
            if '-' not in x:
                pairs.append(x)

        ts_count = 0
        tv_count = 0
        length = len(pairs)

        transitions = ["AG", "GA", "CT", "TC"]
        transversions = ["AC", "CA", "AT", "TA",
                         "GC", "CG", "GT", "TG"]

        for (x, y) in pairs:
            if x + y in transitions:
                ts_count += 1
            elif x + y in transversions:
                tv_count += 1

        p = float(ts_count) / length
        q = float(tv_count) / length
        try:
            d = -0.5 * log((1 - 2 * p - q) * sqrt(1 - 2 * q))
        except ValueError:
            print("Tried to take log of a negative number")
            return None
        return d
