from bin.helper_function import *
from bin.multiAlignGenomeLTRs import RunAndParseClustal
from bin.pathfinder import ProjectStructure
from bin.messager import *
from bin.LTRharvest_loop import LTRharvestRun
from bin.InsertionTimeEstimation import LTR_InsertionTimeCalculator
import argparse
from Bio import SeqIO
import os
from bin import Clustering_LTRs
from bin.screenSelectedLTRbyBLASTgenome import checkLTR_bygBLAST
from bin.selectCluster_for_iRAP import *

class iRAPer():
    def __init__(self, args):
        print("The iRAPer pipeline has been started >>> ")
        self.args = args
        self.genome = args.genome
        self.trna = args.trna
        self.profiles = args.profiles + "/*.hmm"
        self.chunks = self._getChunks()
        self.seq_assign_to_chunks = self._getSeqAssigned()
        self.project_structure = self._getProjectStructure()

    def _getChunks(self):
        showMessageLevel("Fasta file estimation", level=2)
        return genomeSplitterIndex(self.genome)

    def _getProjectStructure(self):
        """
        this woll create all directories needed for the project
        :return:
        """
        showMessageLevel("Project structure generation", level=2)
        project_structure = ProjectStructure(args.out, self.chunks, self.genome)
        print(project_structure)
        return project_structure

    def start_iRAPer(self, chunk, genome_fasta):
        ltr = LTRharvestRun(genome_fasta) #results will be _sorted.gff3 file
        ltr_fasta_path = ltr.runLTRdigest(self.trna, self.profiles) # return path to the 3' LTR fasta generated by LTRdigest
        self.project_structure.addLTRfastaPath(chunk, ltr_fasta_path) # add path to the 3' LTR and 5' LTR


    def createGenomeFastas(self):
        """
        it will create chunks genome fasta in each tmp folder
        :return:
        """
        showMessageLevel("Genome fasta files are creating", level=2)
        chunk_in_process = 0
        genome_fasta = self.project_structure.genomes_tmp[chunk_in_process]
        tmp_fasta = open(genome_fasta, "w")
        files_created = []
        files_created.append(genome_fasta)
        cnt = 1
        for seq in SeqIO.parse(self.genome, "fasta"):
            chunk_current = self.seq_assign_to_chunks[seq.id]  # chunk for this sequence
            if chunk_current != chunk_in_process:
                tmp_fasta.close()
                chunk_in_process = chunk_current
                genome_fasta = self.project_structure.genomes_tmp[chunk_in_process]
                files_created.append(genome_fasta)
                tmp_fasta = open(genome_fasta, "w")
                SeqIO.write(seq, tmp_fasta, "fasta")
                cnt += 1
            else:
                SeqIO.write(seq, tmp_fasta, "fasta")

        showMessageLevel("Number of tmp genome fasta files created {}".format(cnt), level=2)
        return files_created

    def _getSeqAssigned(self):
        """
        assign chunks to indovodual sequences
        :return: dictionary seq.id:chunk
        """
        seq_assigned  = {}
        for i,chunks in enumerate(self.chunks):
            for seq_ids in chunks:
                seq_assigned[seq_ids] = i
        return seq_assigned

    def mergeLTRfasta(self):
        """
        :return: two merged fasta files for 3'LTR and 5'LTR sequences
        """
        ## 3 LTR
        os.system('cat {0} > {1}'.format(' '.join([self.project_structure.ltr_3[i] for i in self.project_structure.ltr_3]),
                  self.project_structure.merged_3))
        ## 5 LTR
        os.system(
            'cat {0} > {1}'.format(' '.join([self.project_structure.ltr_5[i] for i in self.project_structure.ltr_5]),
                                   self.project_structure.merged_5))

    def clustering(self):
        os.system('cdhit -i {0} -o {1} -s 0.4 -aL 0.5 -d 150 -T 0'.format(self.project_structure.merged_3,
                                                                          'cdhit'))
        ## parse the results
        cl = Clustering_LTRs.CdHit_clustering(self.project_structure.cd_hit_results,
                                         self.project_structure.parsed_cd_hit_out)

        showInfoMessage("Maximum sequences ({0}) in cluster {1}".format(len(cl.getMaximumCluster()), cl.getMaximumCluster().id))


    def runBLAST(self):
        os.system("makeblastdb -in {0} -out {0} -dbtype nucl".format(self.project_structure.genome_blastDB))
        #3'-LTR BLAST
        os.system("blastn -query {0} "
                  "-db {1} -outfmt 5 "
                  "-out {2} -evalue 0.000001 "
                  "-window_size 22 -num_threads 20".format(self.project_structure.selection_sequence3_per_cluster,
                                                           self.project_structure.genome_blastDB,
                                                           self.project_structure.BLAST_3LTR_xml))
        #5'-LTR BLAST
        os.system("blastn -query {0} "
                  "-db {1} -outfmt 5 "
                  "-out {2} -evalue 0.000001 "
                  "-window_size 22 -num_threads 20".format(self.project_structure.selection_sequence5_per_cluster,
                                                           self.project_structure.genome_blastDB,
                                                           self.project_structure.BLAST_5LTR_xml))

    def findLTRsInGenome(self):
        checkLTR_bygBLAST(
            self.project_structure.BLAST_3LTR_xml,
            ltr=3)
        checkLTR_bygBLAST(
            self.project_structure.BLAST_3LTR_xml,
            ltr=5)

    def run(self):
        ### step 1: estimate fasta file and calculate number of chunks needed

        ### step 2: create projects structure and corresponding folders and subfolders

        ### step 3:  split genome into fasta file for each chunk
        files_created = self.createGenomeFastas()

        ### step 4: iterate across all files and run LTRharvest and LTRdigest
        for chunk, genome in enumerate(files_created):
            showStep("LTR identification in file: " + genome)
            self.start_iRAPer(chunk, genome)

        ### step 5: merge fasta files for each chunk into 5' LTR and 3' LTR
        showStep("merging of 3' and 5' LTR fasta files from different chunks")
        self.mergeLTRfasta() # results  self.project_structure.merged_3 and self.project_structure.merged_5 fasta files
        showMessageLevel("File {0} and {1} have been created".format(self.project_structure.merged_3,
                                                                     self.project_structure.merged_5),
                         level=1)

        ### step 6: estimate insertion time
        showStep("insertion time estimation")
        # LTR_InsertionTimeCalculator(self.project_structure.merged_3,
        #                             self.project_structure.merged_5,
        #                             self.project_structure.insertion_time_tab)
        showMessageLevel("File {} has been created".format(self.project_structure.insertion_time_tab),
                         level=1)

        ### step 7: clustering of 3’ LTR sequences by CD-HIT
        showStep("LTR clustering")
        self.clustering()

        ### step 8: select clusters and sequences from them
        showStep("Selecting clusters and sequences from them")
        selectClusters_and_LTRs(self.project_structure.parsed_cd_hit_out,
                                self.project_structure.insertion_time_tab,
                                self.project_structure.merged_3,
                                self.project_structure.merged_5,
                                self.project_structure.selection_tab_cluster,
                                self.project_structure.selection_sequence3_per_cluster,
                                self.project_structure.selection_sequence5_per_cluster)
        # it will return fasta of LTRs (single per cluster)  self.selection_sequence_per_cluster = self.root + "/selected_LTR_sequences.fasta"

        ### step 9: BLAST
        showStep("BLAST for selected 3'LTRs and 5'LTRs is going")
        self.runBLAST()

        ### step 10: Collection of LTR similar sequences from genome
        showStep("Collection of LTR similar sequences from whole genome")

        ltr3_collected_sim_seqs = checkLTR_bygBLAST(
            self.project_structure.BLAST_3LTR_xml,
        ltr=3)
        ltr5_collected_sim_seqs = checkLTR_bygBLAST(
            self.project_structure.BLAST_5LTR_xml,
            ltr=5)

        ### step 11: ClustalO multiple alignment of the collected fastas
        showStep('ClustalO multiple alignment of the collected fastas')
        print("\n".join(ltr3_collected_sim_seqs.created_fastas))
        print("////////")
        print("\n".join(ltr5_collected_sim_seqs.created_fastas))
        print("////////")

        for i,ltr3_files in enumerate(ltr3_collected_sim_seqs.created_fastas):
            RunAndParseClustal(ltr3_files, self.project_structure.tmp_folder + "/3LTR_{}.tmp_tab".format(i),ltr=3, run_clustal=True)
            ltr_5_fasta = ltr3_files.replace("::3::","::5::")
            ltr_5_fasta = ltr_5_fasta.replace("LTR3", "LTR5")
            RunAndParseClustal(ltr_5_fasta, self.project_structure.tmp_folder + "/5LTR_{}.tmp_tab".format(i),ltr=5, run_clustal=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Short sample app')

    # Required arguments
    parser.add_argument('genome', help = "Path to the file of genome fasta sequence")
    parser.add_argument('trna', help = "Path to the file with tRNA fasta")
    parser.add_argument('profiles', help = "Path to the folder with .hmm profiles")


    # optional arguments
    parser.add_argument('-o',"--outdir",dest='out',default=".", help = "output directory for iRAPer files")

    #paths to the programmes used by iRAPer
    soft_group = parser.add_argument_group('software paths')
    soft_group.add_argument("-blastn", action="store", default="blastn", help = "path to the blasnt. Default: blastn")

    # version
    parser.add_argument("-v", '--version', action='version',
                        version='1.0')

###########################################################
    args = parser.parse_args()
    iRAPer(args).run()

