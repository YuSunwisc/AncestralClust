### 1_2_initial_sequence_50.py --- This is the script select random initial clust samples(i.e. 50)

# Authors: Yu Sun 
# (Last updated 2021-04-28)
#
# DESCRIPTION: This is the script select random initial clust samples(i.e. 50) and put them into ../data
# folder.

#
# INSTRUCTIONS: Direct to current folder and run
#`python 1_2_initial_sequence_50.py` 
# on your terminal


################################################################################################################################################


### ---Begin script--- ###

## Step 1. Random select use sample package

from Bio import SeqIO
from random import sample
g = open("/Users/shuqi/Desktop/563/Phylo563_FinalProject/data/initail_sequence_50.fasta","w")
with open("/Users/shuqi/Desktop/563/Phylo563_FinalProject/data/sequence_3k.fasta") as f:
    seqs = SeqIO.parse(f, "fasta")
    samps = ((seq.name, seq.seq) for seq in  sample(list(seqs),50))                     # randomly select 50 samples and write into a new fasta file
    for samp in samps:
        g.write(">{}\n{}".format(*samp))
        g.write('\n')
g.close()



### ---End script--- ###

