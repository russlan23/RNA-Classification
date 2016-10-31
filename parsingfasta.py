from intervaltree import Interval, IntervalTree
from collections import namedtuple


###########################
#defining an exon
exon = namedtuple('exon', ['chr_n', 'start', 'end', 'strand'])

def complement(seq):

    complement = { 'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N' }
    complseq = [complement[base] for base in seq]
    return complseq

def reverse_complement(seq):
    seq = list(seq)
    seq.reverse()
    return ''.join(complement(seq))
    
    
###########################
#dealing with bed files
def read_bed(file):
    bedFile = open(file, 'r')
    exons = []

    for line in bedFile:
        chr_n, start, end, strand = line.split(maxsplit=4)
        exons.append(exon(chr_n, int(start), int(end), strand))
    bedFile.close()
    
    return exons

def get_exon_data(exon, seqs):
    data = seqs[exon.chr_n][exon.start:exon.end].upper()

    if exon.strand == '-':
        data = reverse_complement(data)
    return data

def build_interval_trees(exons, seqs):
    trees = {}
    for seq_name, _ in seqs.items():
        trees[seq_name] = IntervalTree()
        
    for exon in exons:
        trees[exon.chr_n][exon.start:exon.end] = exon
    return trees