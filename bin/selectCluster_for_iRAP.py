from collections import defaultdict
from Bio import SeqIO




def selectClusters_and_LTRs(clstr_tab, ins_tab, ltr_fasta3, ltr_fasta5, outFile_tab, outFile_Fasta3, outFile_Fasta5,
                            ltrs_isFull, ltrs_Classification,
                            min_seq_in_cluster = 10,
                            min_percent_young = 70,
                            ins_time_cutoff = 1500000):
    seq_ind3 = SeqIO.index(ltr_fasta3, "fasta")
    seq_ind5 = SeqIO.index(ltr_fasta5, "fasta")

    with open(clstr_tab) as clstr_t, \
        open(ins_tab) as ins_t, \
        open(outFile_Fasta3, "w") as outFasta3, \
            open(outFile_Fasta5, "w") as outFasta5, \
            open(outFile_tab, "w") as out:


        ### estimate age
        ins_time = {} #seq.id = age

        ###classification
        class_te = ltrs_Classification #seq.id = classification

        ###is FUll or not
        isFull_te = ltrs_isFull #seq.id = isFull (True, False)

        for lines in ins_t:
            sp = lines.rstrip().split('\t')
            ins_time[sp[0]] = int(sp[1])
            ids = sp[0].split("_")
            chromosome, start, end = ids[0], ids[-1], ids[-2]

        ### select clusters
        cnt_c_s = 0
        selected = {} # cluster:leading_sequence
        cluster_seqs_age = defaultdict(list) # cluster:[True if age < cutoff or False if > cutoff for each sequence in acluster]
        seq_per_cluster = defaultdict(list) # {cluster:[seq1,seq2,...], }
        for i, lines in enumerate(clstr_t):
            if i != 0:
                sp = lines.rstrip().split('\t')
                seq_id = sp[1]
                num_Seq = sp[-2] # number of sequences in cluster
                if int(num_Seq) >= min_seq_in_cluster:
                    print(lines)
                    cnt_c_s += 1
                    # estimate age of the sequence and True if it < cutoff or False otherwise
                    if ins_time[seq_id] < ins_time_cutoff:
                        cluster_seqs_age[sp[-3]].append(True)
                    else:
                        cluster_seqs_age[sp[-3]].append(False)
                    seq_per_cluster[sp[-3]].append(seq_id)
                    if sp[-1] == "True":
                        selected[sp[-3]] = sp[1]
        print("Number of clusters after filtering by cluster size (>={0}) is {1}".format(min_seq_in_cluster, len(cluster_seqs_age)))

        #######filter by classification and isFull
        TE_class_isFull_selected = {}
        full_per_cluster = {}
        classification_per_cluster = {}
        for cluster in selected:
            isFUll = [isFull_te[id] for id in seq_per_cluster[cluster]]
            full_per_cluster[cluster] = "{0} of {1}".format(isFUll.count('True'), len(isFUll))

            classification = [class_te[id] for id in seq_per_cluster[cluster]]
            classification_per_cluster[cluster] = ";".join([i for i in set(classification)])
            # select cluster if not all TEs in it are truncated and at least on is full

            print(isFUll.count('True'))
            print(classification.count("truncated TE"))
            if isFUll.count('True') > 0 and classification.count("truncated TE") != len(seq_per_cluster[cluster]):
                TE_class_isFull_selected[cluster] = 0

        print("Number of clusters after filtering by TE structure and classification is {0}".format(len(TE_class_isFull_selected)))

        ### select clusters with percentage of you TE < min_percent_young
        cnt = 0
        out.write("\t".join(["Cluster", "Seq.id", "Is.leading", "Cluster size", "Percentage of young TEs", "Percentage of Full TEs", "Classification"]) + "\n")
        for clusters in cluster_seqs_age:
            percent_young = cluster_seqs_age[clusters].count(True)*100/len(cluster_seqs_age[clusters])
            if percent_young > min_percent_young and clusters in TE_class_isFull_selected:
                cnt += 1
                is_leading = True

                # if leading sequence is available for this cluster take one
                if clusters in selected:
                    print("leading sequence: ",selected[clusters])
                    seq_id = selected[clusters]
                # if not take just third one
                else:
                    is_leading = False
                    print("NO leading sequences", seq_per_cluster[clusters][2])
                    seq_id = seq_per_cluster[clusters][2]

                out.write(
                    "\t".join([clusters, seq_id,
                               str(is_leading),
                               str(len(seq_per_cluster[clusters])),
                               str(round(percent_young,2)),
                               full_per_cluster[clusters],
                               classification_per_cluster[clusters]]
                              ) +
                    "\n")
                SeqIO.write(seq_ind3[seq_id], outFasta3, "fasta")
                SeqIO.write(seq_ind5[seq_id], outFasta5, "fasta")
        print("Number of clusters after filtering by cluster size and percentage of young TEs in a cluster is {0}".format(cnt))
#
# selectClusters_and_LTRs(r'C:\Users\Илья\PycharmProjects\iRAPer\cdhit.clstr._parsed.tab',
#                         r'C:\Users\Илья\PycharmProjects\iRAPer\Insertion_time.tab',
#                         r"C:\Users\Илья\PycharmProjects\iRAPer\3_LTR_merged.fasta",
#                         "out",
#                         "outFas.fas")
