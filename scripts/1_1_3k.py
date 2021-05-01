### 1_1_3k.py --- This is the script reduce the whole data set into 3k seqs

# Authors: Yu Sun 
# (Last updated 2021-04-28)
#
# DESCRIPTION: This is the script reduce the whole data set into 3k seqs and put them into ../data
# folder.

#
# INSTRUCTIONS: Direct to current folder and run
#`python 1_1_3k.py` 
# on your terminal


################################################################################################################################################


### ---Begin script--- ###

## Step 1. Random select 3k seqs from whole data set use sample package

from Bio import SeqIO
from random import sample
g = open("/Users/shuqi/Desktop/563/Phylo563_FinalProject/data/sequence_3k.fasta","w")
with open("/Users/shuqi/Desktop/563/Phylo563_FinalProject/data/sequences.fasta") as f:
    seqs = SeqIO.parse(f, "fasta")
    samps = ((seq.name, seq.seq) for seq in  sample(list(seqs),3000))
    for samp in samps:
        g.write(">{}\n{}".format(*samp))
        g.write('\n')
g.close()


### ---End script--- ###


