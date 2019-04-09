from collections import defaultdict
from Bio import SeqIO
import re
class LTR():
    def __init__(self,ID, chromosome, start, end, strand):
        self.ID = ID
        self.start = start
        self.end = end
        self.chromosome = chromosome
        self.strand = strand
        self.outGff3 = ''
        self.features = {"LTRharvest":defaultdict(list), "LTRdigest":defaultdict(list)}

    def addFeature(self, source, name, chromosome, start, end, evalue):
        self.features[source][name].append([chromosome, start, end, evalue])

    def getBasicInfo(self):
        return [self.ID, self.start, self.end]

    def getAllFeatures(self):
        toret = []
        for source in self.features:
            for features in self.features[source]:
                toret.append(features)
        return toret

    def getLTRdigestDomainsCoordinates(self):
        byDomain = []
        for source in self.features:
            for features in self.features[source]:
                for coordinates in self.features[source][features]:
                    domain = features
                    if source != "LTRharvest":
                        domain = features.split("_")[0]

                    coordinates.append(domain)

                    byDomain.append(coordinates)

        return byDomain




    def isFull(self):
        domains_required = ['INT', 'RT', 'GAG', 'RNaseH', "AP"]

        if len(self.features['LTRharvest']['long_terminal_repeat']) != 2 or len(self.features['LTRharvest']["target_site_duplication"]) != 2:
            return False
        #print(self.features['LTRharvest']['long_terminal_repeat'], len(self.features['LTRharvest']['long_terminal_repeat']))

        d = self.splitBydomain()
        for domains in domains_required:
            if domains not in d:
                return False
        return True

    def classify(self, dictionary):
        """
        param dictionary: dictionary from D:\PycharmProjects\MossRepeatomArticle\Scripts\RT_classification_gyDB_hmm.tab
        :
        :return:
        """
        d_source = self.splitBydomain() #domain:source

        all_clas = []
        for rt in d_source["RT"]:
            all_clas.append(dictionary[rt])


        return ",".join([i + ":" + str(all_clas.count(i)) for i in set(all_clas)])

    def getBestHit(self, domain):
        sizes = []
        for features in self.features["LTRdigest"]:
            domain_name = features
            owner_domain = features
            if "_" in features:
                domain_name = features.split("_")[0]
                owner_domain = features.split("_")[1]

            if domain_name == domain:
                #print(self.features["LTRdigest"][features])
                sizes.append([owner_domain,abs(float(self.features["LTRdigest"][features][0][3])), abs(int(self.features["LTRdigest"][features][0][2])) - int(self.features["LTRdigest"][features][0][1])])
        if len(sizes) > 1:
            sizes = [i for i in sizes if i[0] not in ["caulimovirus", "badnavirus", "cavemovirus", "soymovirus"]]
        if len(sizes) != 0:
            min_evalue = min(sizes, key=lambda x: x[1])[1]
            Len = 0
            current_best = ""
            all_best_evalue = []

            for i in sizes:
                if i[1] == min_evalue:
                    all_best_evalue.append(i[0])
                    if i[2] > Len:
                        Len = i[2]
                        current_best = i[0]

            #print(sizes)
            if len(sizes) > 1:
                print(self.ID)
                print(current_best)
                print(":".join([i for i in set(all_best_evalue)]))
            return current_best
        else:
            return "-"


    def splitBydomain(self):
        byDomain = defaultdict(list)
        for source in self.features:
            for features in self.features[source]:
                if source == "LTRharvest":
                    byDomain[features].append(features)
                elif "_" in features:
                    byDomain[features.split("_")[0]].append(features.split("_")[1])
                else:
                    byDomain[features].append(features)
        return byDomain




class LtrDiParser():

    def __init__(self, gff3File, mask_for_chromosome_id = ".[]0"):
        self.gff3File = self.modifyGff3(gff3File)
        self.LTRs = defaultdict(LTR)
        self.mask_for_chromosome_id = mask_for_chromosome_id ##where split and which index to extract chromosome id from sequence id
        self.run()

    def run(self):
        self.__readGff()

    def __getChromosomeId(self, seqId):
        return seqId.split(self.mask_for_chromosome_id.split('[]')[0])[int(self.mask_for_chromosome_id.split('[]')[1])]

    def __readGff(self):
        with open(self.gff3File) as gff_file:
            start_new_LTR = False
            for num,lines in enumerate(gff_file):
                sp = lines.rstrip().split("\t")
                if not lines.startswith("##gff-version"):
                    if lines.startswith("###") or num == 1:
                        start_new_LTR = True
                    elif start_new_LTR:
                        start_new_LTR = False
                        newLTR = LTR(sp[-1].split("=")[1], sp[0].split(' ')[0], sp[3], sp[4], sp[6])
                        self.LTRs[newLTR.ID] = newLTR
                    else:
                        if sp[1] == "LTRharvest" or ";name=" not in sp[-1]:
                            featureName = sp[2]
                        else:
                            featureName = sp[-1].split(";name=")[1]
                        self.LTRs[newLTR.ID].addFeature(sp[1],featureName, sp[0], sp[3], sp[4], sp[5])

        print("Number of LTRs found: ", len(self.LTRs))


    def getAllfeatureNames(self):
        fn = []
        domain = []
        domain_presents = []
        for ltrs in self.LTRs:
            fn += self.LTRs[ltrs].getAllFeatures()
            dLTRs = self.LTRs[ltrs].splitBydomain()
            domain += dLTRs
            domain_presents.append(",".join(sorted(dLTRs)))

        # print(set(fn))
        # print(len(set(fn)))
        # print(set(domain))
        #
        # print(len(set(domain_presents)))
        # for dP in set(domain_presents):
        #     print(dP, domain_presents.count(dP))

    def getAllFull(self):
        cnt = 0
        for ltrs in self.LTRs:
            if self.LTRs[ltrs].isFull():
                cnt += 1
        print(cnt)

    def getFastaFullLtrs(self, genomeFasta):
        cnt_t = 0
        per_chromosome = defaultdict(list)
        for ltrs in self.LTRs:
            if self.LTRs[ltrs].isFull():
                per_chromosome[self.LTRs[ltrs].chromosome].append(self.LTRs[ltrs].getBasicInfo())
                cnt_t += 1

        cnt = 0

        with open("Full_LTRs_transposons_{}".format(genomeFasta.split("/")[-1]), "w") as outFile:
            for seq in SeqIO.parse(genomeFasta, "fasta"):
                if seq.id in per_chromosome:
                    chromosome_seq = str(seq.seq)
                    if seq.id == "scaffold_49":
                        print(chromosome_seq)
                    for ltrs in per_chromosome[seq.id]:
                        #print(ltrs)
                        start, end = min([int(ltrs[1]), int(ltrs[2])]), max(([int(ltrs[1]), int(ltrs[2])]))
                        outFile.write(">" + ltrs[0] + " {2}:{0}_{1}".format(start-1,end, seq.id) + "\n" +  chromosome_seq[start-1:end] + "\n")
                        cnt += 1

        print("{0} of {1} sequences in the file".format(cnt, cnt_t))


    def getClassification(self):
        with open(r"RT_classification_gyDB_hmm.tab") as classTab, \
            open("fullLTRclassification{}".format(self.gff3File.split("/")[-1]), "w") as outfile:
            class_d = {}
            class_d['micropia'] = "Ty3/Gypsy"
            class_d['v'] = "Ty3/Gypsy"
            class_d['412'] = "Ty3/Gypsy"
            class_d['17'] = "Ty3/Gypsy"
            class_d['a'] = "Ty3/Gypsy"
            class_d['b'] = "Ty3/Gypsy"
            class_d['codi'] = "Ty1/Copia"
            class_d["galea"] = "Ty1/Copia"


            for lines in classTab:
                sp = lines.rstrip().split()
                class_d[sp[0]] = sp[1]
            for ltrs in self.LTRs:
                base = "\t".join([ltrs,
                                  self.__getChromosomeId(self.gff3File.split("/")[-1].replace('_',"|")),
                                  self.LTRs[ltrs].start, self.LTRs[ltrs].end ])
                if  self.LTRs[ltrs].isFull():
                    print(self.LTRs[ltrs].ID)
                    outfile.write(base + "\t" + self.LTRs[ltrs].classify(class_d) + "\t" + class_d[self.LTRs[ltrs].getBestHit("RT")] + "\n")
                else:
                    outfile.write(base + "\t" + 'truncated TE' + "\t" + 'truncated TE' + "\n")

    def getBEDfileDomains(self, from0 = True):
        with open("{}.bed".format(self.gff3File[:-4]), "w") as outfile:
            for re in self.LTRs:
                per_domain = defaultdict(list)
                if self.LTRs[re].isFull():
                    al_domains_matrix = self.LTRs[re].getLTRdigestDomainsCoordinates()

                    for al_domains in al_domains_matrix:
                        per_domain[al_domains[-1]].append(al_domains)

                    #print(per_domain["LTR_retrotransposon"][0])
                    start = int(per_domain["LTR_retrotransposon"][0][1])
                    for domains_unqiue in per_domain:
                        if domains_unqiue != "long_terminal_repeat":
                            longest_hit = max(per_domain[domains_unqiue], key=lambda x:abs(int(x[1]) - int(x[2])))
                            if from0:
                                longest_hit[1] = str(int(longest_hit[1]) - start)
                                longest_hit[2] = str(int(longest_hit[2]) - start)
                            outfile.write(re + "\t" + "\t".join(longest_hit) + "\n")
                        else:
                            for ltrs in per_domain[domains_unqiue]:
                                if from0:
                                    ltrs[1] = str(int(ltrs[1]) - start)
                                    ltrs[2] = str(int(ltrs[2]) - start)
                                outfile.write(re + "\t" + "\t".join(ltrs) + "\n")


    def gff3Tobed(self):
        with open("{}.bed".format(self.gff3File[:-4]), "w") as outfile:
            for re in self.LTRs:
                per_domain = defaultdict(list)
                al_domains_matrix = self.LTRs[re].getLTRdigestDomainsCoordinates()

                for al_domains in al_domains_matrix:
                    per_domain[al_domains[-1]].append(al_domains)

                #print(per_domain["LTR_retrotransposon"][0])
                start = int(per_domain["LTR_retrotransposon"][0][1])
                for domains_unqiue in per_domain:
                    if domains_unqiue != "long_terminal_repeat":
                        longest_hit = max(per_domain[domains_unqiue], key=lambda x:abs(int(x[1]) - int(x[2])))
                        outfile.write("\t".join(longest_hit) + "\t" + re + "\n")
                    else:
                        for ltrs in per_domain[domains_unqiue]:
                            outfile.write("\t".join(ltrs) + "\t" + re + "\n")



    def modifyGff3(self, LTRharvest_output):
        """
        instead of real chromosome name LTharvest put Seq.. names. This function change chromosome and scaffolds names
        :param LTRharvest_output:
        :return:
        """
        l1_new_names = [] # seq...
        l2_real_names = [] # Chr...
        d = {} # number in gff3 annotation : real chromosome name
        real_name_patter = re.compile("^#[^#]")
        outfile_id = LTRharvest_output + '_modified'
        real_names_start = False
        with open(LTRharvest_output) as infile, open(outfile_id, "w") as outfile:
            for num, lines in enumerate(infile):
                if num != 0:
                    if lines.startswith("##s"):
                        l1_new_names.append(lines.split(" ")[3].split("seq")[-1])
                    elif real_name_patter.search(lines):
                        l2_real_names.append(lines.rstrip()[1:])
                    elif lines.startswith("###"):
                        outfile.write(lines)
                    else:
                        if not real_names_start:
                            l1_new_names = [int(i) for i in l1_new_names]
                            for i,nam in enumerate(sorted(l1_new_names)):
                                d["seq" + str(nam)] =  l2_real_names[i]
                            real_names_start = True

                        sp = lines.split("\t")
                        sp[0] = d[sp[0]]
                        new_line = "\t".join(sp)
                        outfile.write(new_line)
                else:
                    outfile.write("###")


        return outfile_id
    def changeIDseqs(self, fastafile):
        with open("modifiedIDs_" + fastafile, "w") as outfile:
            for seq in SeqIO.parse(fastafile, "fasta"):
                seq.id = seq.id + "_" + "_".join(seq.description.split(" "))

                SeqIO.write(seq, outfile, "fasta")

    def findOverlap(self, chromosome, start, end):
        for ltrs in self.LTRs:
            ltrs = self.LTRs[ltrs]
            ltrs.start, ltrs.end = int(ltrs.start), int(ltrs.end)
            #print(ltrs.ID, ltrs.chromosome)
            if ltrs.chromosome == chromosome:
                if not(start > ltrs.end and end > ltrs.end) and not (start < ltrs.start and end < ltrs.start):
                    print(ltrs.ID)
    #
# LD = LtrDiParser(r"C:\Users\Илья\PycharmProjects\iRAPer\genome_chunk_1.fasta.idx_LtrDi.gff3")
# LD.findOverlap('CP027625.1', 9055592, 9060554)
# #LD.gff3Tobed()
# LD.getBEDfileDomains()
# #LD.getAllfeatureNames()
#LD.getClassification()
# # #LD.getFastaFullLtrs("/home/ilia/RNA_seq_reads_moss/N/Ppatens_318_v3.fa")