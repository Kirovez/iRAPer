from Bio import SeqIO
from bin.messager import *
def genomeSplitterIndex(fasta, size = 100):
    """
    function to estimate number of chunks genome split to ..
    :param fasta: genome file
    : size : size of a chunk in Mb
    :return: list of list with sequence ids for each chuck
    """
    cutoff = size * 1000000
    chunk_list = [[]]
    accumulated_size = 0
    for seq in SeqIO.parse(fasta, "fasta"):
        accumulated_size += len(seq.seq)
        chunk_list[-1].append(seq.id)
        if accumulated_size > cutoff:
            accumulated_size = 0
            chunk_list.append([])
        else:
            chunk_list[-1].append(seq.id)

    if not chunk_list[-1]:
        chunk_list = chunk_list[:-1]

    if len(chunk_list) > 1:
        showInfoMessage("Genome file was splitted into {} chunks...".format(len(chunk_list)))
    return chunk_list

