import os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from iRAPer import primerTest


class parseMafft():
    def __init__(self, fasta, ltr = 3, run_mafft=True):
        self.fasta_file = fasta
        self.mafft_fasta = '{0}.cls_fasta'.format(self.fasta_file)
        self.score = {} # {0: {'a': 107, 't': 0, 'g': 0, 'c': 0, '-': 12}, 1: {'a': 0, 't': 108, 'g': 0, 'c': 0, '-': 11}, 2: {'a': 0, 't': 0, 'g': 0, 'c': 112, '-': 7},
        self.seq_length = 0
        self.sequence_list = []
        self.sequence_count = 0
        self.run_mafft = run_mafft
        self.runMafft()
        self.window = 20
        self.outTable = open(self.mafft_fasta + "oligo_primers.tab", 'w')
        self.writeHeader()
        self.ltr = ltr
        self.main()

    def runMafft(self):
        if self.run_mafft:
            os.system('clustalo --threads 20 -i {0} > {1}'.format(self.fasta_file, self.mafft_fasta))


    def main(self):
        self.score = self.count()
        con_bases, cons_score = self.build_consensus(self.score)
        self.find_regions_for_primers(con_bases, cons_score)

    def count(self):
        score = {}
        for i, seq in enumerate(SeqIO.parse(self.mafft_fasta, 'fasta')):
            self.sequence_list.append(str(seq.seq))
            self.sequence_count += 1
            for n_base, bases in enumerate(str(seq.seq)):
                if i == 0:
                    self.seq_length = len(seq.seq)
                    score[n_base] = {'A':0, 'T':0, 'G':0, 'C':0, '-':0, 'N':0}
                score[n_base][bases] += 1
        #print(score)

        return score

    def build_consensus(self, score_dic):
        cons_score = []
        cons_bases = []
        for i in range(self.seq_length):
            base = max(score_dic[i], key=score_dic[i].get)
            score_base = score_dic[i][base]
            cons_score.append(score_base)
            cons_bases.append(base)
        return [cons_bases, cons_score]

    def absolutFrequency(self, primer):
        """
        It will count absolute frequency of the primers in sequences used for alignment
        :param primer: primer_sequence
        :return: frequency
        """
        cnt = 0
        for seq in self.sequence_list:
            if primer in str(seq):
                cnt += 1

        return [round(cnt/self.sequence_count,2), '{0}/{1}'.format(cnt, self.sequence_count)]


    def writeHeader(self):
        self.outTable.write('\t'.join(['File',
                                       'Start',
                                       'End',
                                       'Min score',
                                      'Mean score',
                                       'First base score',
                                       'Last base score',
                                       'Sequence',
                                       'Primer',
                                       'Frequence of the sequence',
                                       'Absolute number of seqs with the primer',
                                       "GC",
                                       "Tm",
                                       'LTR']) + "\n")

    def find_regions_for_primers(self, cons_base, cons_score):
        max_score = self.sequence_count*self.window
        cutoff = 0.5

        for start, score_value in enumerate(range(len(cons_score) - self.window)):
            current_window = cons_score[start:start + self.window]
            sequence = "".join(cons_base[start:start + self.window])

            ##scores
            min_score_in_window = round(min(current_window)/self.sequence_count,1)
            mean_score_in_window = round((sum(current_window) / self.window)/self.sequence_count)
            first_base_score = round(current_window[0]/self.sequence_count,1)
            last_base_score = round(current_window[-1] / self.sequence_count, 1)
            if min_score_in_window > cutoff and sequence.count("-") == 0:
                # print(current_window)
                # print(sequence)
                # print("start:{0}\nend:{1}\n".format(start + 1, start + self.window))
                primer = sequence
                if self.ltr == 5: # write reverse comlement sequence as a primer
                    primer = str(Seq(primer).reverse_complement())
                Tm = round(primerTest.calculateTm(primer), 2)
                if (primer[-1] == 'G' or primer[-1] == 'C') and Tm > 53:
                    self.outTable.write('\t'.join([str(i) for i in [self.mafft_fasta,
                                                   start + 1,
                                                   start + self.window,
                                                   min_score_in_window,
                                                   mean_score_in_window,
                                                   first_base_score,
                                                   last_base_score,
                                                   sequence,
                                                   primer,
                                                   *self.absolutFrequency(sequence),
                                                                    GC(Seq(primer)),
                                                                    Tm,
                                                   str(self.ltr) + "LTR"]]) + "\n")


#parseMafft(r'5LTR_00_unplaced_join_73N_24795102_24803627Cluster25010.fasta', run_mafft=False)