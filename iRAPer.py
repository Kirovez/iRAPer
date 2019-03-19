from bin.helper_function import *
from bin.pathfinder import ProjectStructure
from bin.messager import *
from bin.LTRharvest_loop import LTRharvestRun
import argparse
from Bio import SeqIO

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
        showMessageLevel("Project structure generation", level=2)
        project_structure = ProjectStructure(args.out, self.chunks)
        print(project_structure)
        return project_structure

    def start_iRAPer(self, genome_fasta):
        ltr = LTRharvestRun(genome_fasta) #results will be _sorted.gff3 file
        ltr.runLTRdigest(self.trna, self.profiles)

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

    def run(self):
        files_created = self.createGenomeFastas()
        # iterate across all files and run LTRharvest
        for genome in files_created:
            showStep("LTR identification in file: " + genome)
            self.start_iRAPer(genome)

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

