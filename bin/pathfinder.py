import os
"""
functions to make folder structure for iRAPer output
"""

class ProjectStructure:
    def __init__(self,outDir, chunks):
        """
        :param outDir: path to out directory
        :param chunks: list od seq.id chunks returned by genomeSplitterIndex function
        :return: dictionary
        """
        self.root = self._isExistCreate(outDir)
        self.chunks = chunks
        self.LTRh_tmp = self._getLTRharvestOutDir() # all LTRharvest files dictionary where chunk is a key
        self.genomes_tmp = self._getGenomeFile() #ids for all genome fasta where chunk is a key\
        #self.LTRharvest_out_files_suffixes = ['_LTRs_inner.fas','_LTRs.gff3','.fa']
        self.LTRharvest_gff3 = '_LTRs.gff3'
    def _getLTRharvestOutDir(self):
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


    def __str__(self):
        return "root directoty: {0} \n" \
               "tmp files \n{1}" \
               "Table with primers \n".format(self.root,
                         "\n".join(["\t\t\t" + self.LTRh_tmp[i] for i in self.LTRh_tmp])
                         )
