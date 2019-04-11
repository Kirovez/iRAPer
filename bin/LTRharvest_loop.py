import os
from bin.messager import *
from Bio import SeqIO
from bin.LtrDiParser import *

class LTRharvestRun():
    def __init__(self,genome_fasta,skip=False):
        self.genome_fasta = genome_fasta
        self.index_name_root = genome_fasta + ".idx"
        self.index_id = self.index_name_root.split('/')[-1]
        self.gff3 = self.index_name_root + '_LTRs.gff3'
        self.skip = skip
        self.runLTRharvest()
        self.sorted_gff3 = self.sortGff3(self.gff3)

    def isSingleSeqInFile(self):
        """
        When only one sequence in fasta file then the gff3 generated will be without real name in header
        Therefore this nreal name of a sequence must be provided for LtrDiParser as sequence_name
        :return:
        """
        cnt = 0
        id = ''
        for seq in SeqIO.parse(self.genome_fasta, 'fasta'):
            cnt += 1
            id = seq.id

        if cnt == 1:
            return id
        else:
            return False

    def runLTRharvest(self):
        ## three outfiles differ in the format wrtiting: .fa, .gff3 and .fas
        print('gt suffixerator -db {0} -indexname {1} -tis -suf -lcp -des -ssp -sds -dna'.format(self.genome_fasta, self.index_name_root ))

        if not self.skip:
            os.system('gt suffixerator -db {0} -indexname {1} -tis -suf -lcp -des -ssp -sds -dna'.format(self.genome_fasta, self.index_name_root ))
            os.system('gt ltrharvest -index {0} -v  -out {0}_LTRs.fa -gff3 {0}_LTRs.gff3 -outinner {0}_LTRs_inner.fas'.format(self.index_name_root))

    def sortGff3(self,gff3):
        ## sort gff3 ###
        showMessageLevel("Gff3 sorting" ,level=1)
        if not self.skip:
            os.system('gt gff3 -sort {0} > {0}_sorted.gff3'.format(gff3))
        return self.gff3 + '_sorted.gff3'

    def getLTR_from_LTRharvesGff3(self):
        """
        :return: will create fasta files for 3' and 5' LTRs and return the full path to them
        """
        single_id = self.isSingleSeqInFile()
        if single_id:
            LD = LtrDiParser(self.sorted_gff3, sequence_name=single_id).get_LTRs_fasta(self.genome_fasta)
        else:
            LD = LtrDiParser(self.sorted_gff3).get_LTRs_fasta(self.genome_fasta)

        return LD

    def runLTRdigest(self, trna_fasta, domains_folder):
        index_genome = self.index_name_root
        sorted_gff3 = self.sorted_gff3
        showStep("LTRdigest is running")
        ltrDigest_command = 'gt ltrdigest -force -outfileprefix {0} -o {0}_LtrDi.gff3 -pptlen 10 30 ' \
                            '-pbsoffset 0 3 -trnas {1} -hmms {2} -pdomevalcutoff 0.001 {3} {4}'.format(
            self.index_name_root,
            trna_fasta,
            domains_folder,
            sorted_gff3,
            index_genome
        )
        print(ltrDigest_command)
        if not self.skip:
            os.system(ltrDigest_command)

        single_id = self.isSingleSeqInFile()

        if single_id:
            LD = LtrDiParser('{0}_LtrDi.gff3'.format(self.index_name_root), sequence_name=single_id)
        else:
            LD = LtrDiParser('{0}_LtrDi.gff3'.format(self.index_name_root))
        classificati_file = LD.getClassification()
        return classificati_file


