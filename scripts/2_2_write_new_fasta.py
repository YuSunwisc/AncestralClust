### 2_2_write_new_fasta.py --- This script writes pairwise fasta data with all centroids

# Authors: Yu Sun 
# (Last updated 2021-04-29)
#
# DESCRIPTION: This script writes pairwise fasta data with all centroids.
#
# INSTRUCTIONS: Run in py3.


################################################################################################################################################


### ---Begin script--- ###

## Step 1. Write centriod sequence into centroid.fasta file

for i in range(1,6):
    allseq = "/Users/shuqi/Desktop/563/Phylo563_FinalProject/data/sequence_3k.fasta"
    #allseq = "/Users/shuqi/Desktop/563/Phylo563_FinalProject/data/example_input.fasta"

    
    centroid_fasta = open("/Users/shuqi/Desktop/563/Phylo563_FinalProject/data/initial_clust_"+str(i)+"_centroid.fasta","w")  #centroid.fasta
    centroid_txt = open("/Users/shuqi/Desktop/563/Phylo563_FinalProject/data/initial_clust_"+str(i)+"_centroid.txt","r")       #centroid.txt file
    centroid_name = centroid_txt.readline().rstrip('\n')            # get centroid name
    with open(allseq,"r") as input:
        iter_input = iter(input)                                    # convert input to interable type so can use next for lines
        for idx,line in enumerate(iter_input):                      # enumerate can help get index and value at same type, but unnecessary here
            if centroid_name in line:                               # find out the name
                centroid_fasta.write(line)
                next_line = next(iter_input)
                while not ">" in next_line:                         # copy all the sequence until see next ">" as header
                    centroid_fasta.write(next_line)
                    next_line = next(iter_input)

    centroid_fasta.close()
    centroid_txt.close()
    


## Step 2. Write centriod sequence & new sample sequence into new_sample.fasta file

from Bio import SeqIO
from random import sample

    
with open("/Users/shuqi/Desktop/563/Phylo563_FinalProject/data/sequence_3k.fasta") as f: # all seq file
    seqs = SeqIO.parse(f, "fasta")
    samps = ((seq.name, seq.seq) for seq in  sample(list(seqs),3000))
    idx = 1
    for samp in samps:
        for i in range(1,6):
            new_fasta = open("/Users/shuqi/Desktop/563/Phylo563_FinalProject/data/pairwise with centroids/"+str(i)+"/new_sample_"+str(idx)+".fasta","w")  #new_sample.fasta
            new_fasta.write(">{}\n{}".format(*samp))                                                                     # write new sample seq into fasta
            new_fasta.write('\n')
            centroid_fasta = open("/Users/shuqi/Desktop/563/Phylo563_FinalProject/data/initial_clust_"+str(i)+"_centroid.fasta","r")  #centroid.fasta
            centroid_seq = centroid_fasta.readlines()                                                                       #write centroid seq into fasta
            for line in centroid_seq:
                new_fasta.write(line)
        idx = idx+1
            




### ---End script--- ###







#######################################################  Appendix  ####################################################################################


# Example of xxx
##############################################


##############################################
