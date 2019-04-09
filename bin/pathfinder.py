import os
"""
functions to make folder structure for iRAPer output
"""

class ProjectStructure:
    def __init__(self,outDir, chunks, fasta):
        """
        :param outDir: path to out directory
        :param chunks: list od seq.id chunks returned by genomeSplitterIndex function
        :return: dictionary
        """
        self.root = self._isExistCreate(outDir)
        self.outTable = self.root + "/iRAPer_results.tab"
        self.chunks = chunks
        self.tmp_folder = self._isExistCreate(outDir + "/tmp")
        self.LtrDi_files = []
        self.LTRh_tmp = self._getLTRharvestOutDir() # all LTRharvest files dictionary where chunk is a key
        self.genomes_tmp = self._getGenomeFile() #ids for all genome fasta where chunk is a key\
        self.LTRharvest_gff3 = '_LTRs.gff3'
        self.ltr_3 = {} # dictionary created after LTRdigest will be run. {chunk:ltr_3 fasta file, ... }
        self.ltr_5 = {} # dictionary created after LTRdigest will be run. {chunk:ltr_5 fasta file, ... }
        self.merged_3 = self.root + "/3_LTR_merged.fasta"
        self.merged_5 = self.root + "/5_LTR_merged.fasta"
        self.insertion_time_tab = self.root + "/Insertion_time.tab"
        self.cd_hit_results = self.root + "/cdhit.clstr"
        self.parsed_cd_hit_out = self.cd_hit_results + '._parsed.tab'
        self.selection_tab_cluster = self.root + "/selected_clusters.info"
        self.selection_sequence3_per_cluster = self.root + "/selected_LTR_3_sequences.fasta"
        self.selection_sequence5_per_cluster = self.root + "/selected_LTR_5_sequences.fasta"
        self.genome_blastDB = fasta
        self.outBLASTdir = self._isExistCreate(self.root + "/BLAST_out")
        self.BLAST_3LTR_xml = self.outBLASTdir  + "/LTR3_vs_Genome_blast5fmt.xml"
        self.BLAST_5LTR_xml = self.outBLASTdir + "/LTR5_vs_Genome_blast5fmt.xml"


    def _getLTRharvestOutDir(self):
        """
        :return: dictionary where key is chunk number and
        value is full path to the LTRharvest output folder for this chunk.
        like: root/chunk/tmp_LTRharvest
        """
        name_getLTRharvestOutDir = "tmp_LTRharvest"
        toret = {}
        if len(self.chunks) == 1:
            toret =  {0:self.root + "/tmp_LTRharvest"}
        else:
            for i, chunks in enumerate(self.chunks):
                self._isExistCreate(self.root + "/" + str(i))
                toret[i] = self.root + "/" + str(i) + "/{}".format(name_getLTRharvestOutDir)

        for tmp in toret:
            self._isExistCreate(toret[tmp])
        return toret

    def _getGenomeFile(self):
        toret = {}
        if len(self.chunks) == 1:
            toret =  {0:self.LTRh_tmp[0] + "/genome.fasta"}
        else:
            for i, chunks in enumerate(self.chunks):
                toret[i] = self.LTRh_tmp[i]  + "/genome_chunk_{}.fasta".format(i)
        return toret

    def _isExistCreate(self, directory):
        if not os.path.isdir(directory):
            os.system('mkdir {}'.format(directory))
        return directory

    def addLTRfastaPath(self, chunk, ltr_fasta):
        self.ltr_3[chunk] = self.genomes_tmp[chunk] + ".idx_3ltr.fas"
        self.ltr_5[chunk] = self.genomes_tmp[chunk] + ".idx_5ltr.fas"


    def __str__(self):
        return "root directoty: {0} \n" \
               "tmp files \n{1}\n" \
               "Table with primers \n".format(self.root,
                         "\n".join(["\t\t\t" + self.LTRh_tmp[i] for i in self.LTRh_tmp])
                         )
