from intervaltree import Interval, IntervalTree

##############################
#Gloabl Variables:
treeCircRNAExons= IntervalTree()

##############################


def parseExonsIntoTree(file, strand, chrom):
    
    print("Parsing allExons into Interval Tree")
    bedFile = open(file, 'r')
    tree= IntervalTree()
    for line in bedFile:
        chr_n, start, end, strand = line.split(maxsplit=4)
        if chr_n ==chrom and strand == strand: 
            start, end = int(start), int(end)
            if end >start:
                 tree.add(Interval(start,end))
    bedFile.close()
    return tree

##############################

def parseCircRNAExonsIntoTree(file, strand, chrom):
    
    print("Parsing exons frpom circRNA into sets")
    bedFile = open('data/hsa_hg19_Rybak2015.bed', 'r')
    
    startsCircExons=[]
    endsCircExons=[]
    
    for line in bedFile:
        chr_n, start, end, gene_name, _, strand, _, _, _, _, sizes, starts = line.split(maxsplit=12)
        start, end = int(start), int(end)
        starts = [int(n) for n in starts.split(',')]
        sizes = [int(n) for n in sizes.split(',')]
        
        if chr_n == 'chr1' and strand == '+':
            for i in range(0, len(starts) ):
              if start + starts[i] not in startsCircExons:
                  startsCircExons.append(start + starts[i])
                  endsCircExons.append(start + starts[i] + sizes[i])
               
                                                                                                          
    bedFile.close()


##############################
allExons='data/all_exons.bed'
positiveChr1Tree=parseExonsIntoTree(allExons,"+","chr1")
parseCircRNAExonsIntoTree

##############################
print("Parsing exons frpom circRNA into sets")
bedFile = open('data/hsa_hg19_Rybak2015.bed', 'r')

startsCircExons=[]
endsCircExons=[]

for line in bedFile:
    chr_n, start, end, gene_name, _, strand, _, _, _, _, sizes, starts = line.split(maxsplit=12)
    start, end = int(start), int(end)
    starts = [int(n) for n in starts.split(',')]
    sizes = [int(n) for n in sizes.split(',')]
    
    if chr_n == 'chr1' and strand == '+':
        for i in range(0, len(starts) ):
          if start + starts[i] not in startsCircExons:
              startsCircExons.append(start + starts[i])
              endsCircExons.append(start + starts[i] + sizes[i])
           
                                                                                                      
bedFile.close()


#############################################
for x in range(len(startsCircExons)):
    positiveChr1Tree.remove_overlap(startsCircExons[x], endsCircExons[x])
      
##############################
#creating negative examples 
newbed=open('data/negative_chr1_plus.bed','w')

for interval_obj in positiveChr1Tree:
   set=positiveChr1Tree[interval_obj]
   for x in set:
       newbed.write('{0}\t{1:d}\t{2:d}\t{3:s}\n'.format('chr1', x.begin, x.end,"+"))
newbed.close()

###############################

#creating negative examples 
newbed=open('data/positif_chr1_plus.bed','w')

for x in range(len(startsCircExons)):
    newbed.write('{0}\t{1:d}\t{2:d}\t{3:s}\n'.format('chr1', startsCircExons[x], endsCircExons[x], "+"))
newbed.close()