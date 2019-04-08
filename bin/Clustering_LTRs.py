class Cluster():
    def __init__(self, id):
        self.id = id
        self.sequences_ids_in_cl = {} ## seq_id : [chromosome_id, start, end]

    def __len__(self):
        return len(self.sequences_ids_in_cl)

    def isLeadingInCluster(self, seq_id):
        if "*" in seq_id:
            return True
        else:
            return False

    def append_sequence(self, seq_id):
        id = seq_id.split("...")[0].split('.')[0]
        start, end = seq_id.split('_')[1], seq_id.split('_')[2].split("...")[0]
        self.sequences_ids_in_cl[id] = [seq_id.split("...")[0], start, end, self.isLeadingInCluster(seq_id)]

    def getLines(self):
        cl = []
        for sequences in self.sequences_ids_in_cl:
            cl.append("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}".format(sequences,
                                               self.sequences_ids_in_cl[sequences][0],
                                               self.sequences_ids_in_cl[sequences][1],
                                               self.sequences_ids_in_cl[sequences][2],
                                               self.id,
                                               str(len(self.sequences_ids_in_cl)),
                                                self.sequences_ids_in_cl[sequences][3]))
        if len(cl) > 1:
            cl = "\n".join(cl)
        else:
            cl = cl[0]

        return cl + "\n"

class CdHit_clustering():
    def __init__(self, CLSTR_FILE,  outFile):
        self.CLSTR_FILE = CLSTR_FILE
        self.Clusters_list = []
        self.outTab = outFile
        self.parseResults()
        self.getOutFile()

    def addToLastCluster(self, sequence_id):
        self.Clusters_list[-1].append_sequence(sequence_id)

    def __len__(self):
        return len(self.Clusters_list)

    def getMaximumCluster(self):
        return max(self.Clusters_list, key= lambda x:len(x))

    def parseResults(self):
        with open(self.CLSTR_FILE) as infile:
            for lines in infile:
                if lines.startswith('>'):
                    cl_id = lines.rstrip()[1:].split(' ')[-1]
                    self.Clusters_list.append(Cluster(cl_id))
                else:
                    self.addToLastCluster(lines.rstrip().split(">")[1])

    def getOutFile(self):
        with open(self.outTab, 'w') as outFile:
            outFile.write("Chromosome\tSequence_id\tStart\tEnd\tCluster\tNum.Seq.In.Cluster\tIs.lead.In.Cluster\n")
            for clusters in self.Clusters_list:
                outFile.write(clusters.getLines())

#cl = CdHit_clustering('3_LTR_merged.clstr')
#print("Maximum sequences ({}) in cluster ".format(len(cl.getMaximumCluster())), cl.getMaximumCluster().id)