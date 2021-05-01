### 0_preparation.sh --- This is the terminal file for all the preparation

# Authors: Yu Sun 
# (Last updated 2021-04-27)
#
# DESCRIPTION: This script includes all the prep steps for data and data analysis
#
# INSTRUCTIONS: Direct to current folder
# `bash 0_preparation.sh $1`
# $1 is the desired initial cluster #


################################################################################################################################################


### ---Begin script--- ###

## Step 1. Install pip3 and Biopy
brew install python3
pip3 install biopython



## Step 2. Prep all empty fasta file

touch ../data/sequence_3k.fasta             # sub samples from the whole dataset
touch ../data/initail_sequence_50.fasta     # randomly select ~1.5% from the sub samples according to the algorithm
touch ../data/pairwisc.fasta                # for writing pairwise seq and calculate distance in R

mkdir ../data/"pairwise with centroids"     # make directory and new fasta file for new sample and centroids
for i in {1.."$1"}                          
do
    mkdir ../data/"pairwise with centroids"/"$i"
    for j in {1..3000}
    do
        touch ../data/"pairwise with centroids"/"$i"/new_sample_"$j".fasta
    done 
done


## Step 3. Prep all empty txt file for initial cluster data

mkdir ../data/"final_clust"     # make directory and final clust data
for i in {1.."$1"}
do
    touch ../data/initial_clust_"$i".txt  # create txt file for each clust centroid header
    touch ../data/initial_clust_"$i"_centroid.txt  # create txt file for each clust centroid header
    touch ../data/initial_clust_"$i"_centroid.fasta  # create fasta file for each clust centroid seq
    touch ../data/initial_clust_"$i"_centroid_pls_new.fasta  # create fasta file for each clust centroid seq and new samples
    touch ../data/"final_clust"/final_clust_"$i".txt  # create txt file for final class member without duplication

done

### ---End script--- ###







#######################################################  Appendix  ####################################################################################


# Example of xxx
##############################################


##############################################


